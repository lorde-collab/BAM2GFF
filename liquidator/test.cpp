#include "gtest/gtest.h"

#include "score_matrix.h"
#include "detail/score_matrix_detail.h"

using namespace liquidator;

const std::array<double, AlphabetSize> uniform_bg = {.25, .25, .25, .25};

TEST(ScoreMatrix, read_pwm_matrix)
{
    const std::string input_str = R"(MEME version 4

ALPHABET= ACGT

strands: +

Background letter frequencies
A 0.29 C 0.21 G 0.21 T 0.29

MOTIF JASPAR2014.MA0107.1 RELA

letter-probability matrix: alength= 4 w= 10 nsites= 18 E= 0
  0.000000        0.222222        0.611111        0.166667
  0.000000        0.000000        0.944444        0.055556
  0.000000        0.000000        1.000000        0.000000
  0.611111        0.000000        0.388889        0.000000
  0.555556        0.166667        0.222222        0.055556
  0.111111        0.000000        0.000000        0.888889
  0.000000        0.000000        0.000000        1.000000
  0.000000        0.111111        0.000000        0.888889
  0.000000        1.000000        0.000000        0.000000
  0.000000        1.000000        0.000000        0.000000)";

    std::istringstream ss(input_str);
    std::vector<detail::PWM> pwms = detail::read_pwm(ss);
    ASSERT_EQ(1, pwms.size());
    
    const auto& name = pwms[0].name;
    EXPECT_EQ("JASPAR2014.MA0107.1", name);

    const auto& matrix = pwms[0].matrix;
    ASSERT_EQ(10, matrix.size());
    const auto& row = matrix[0];
    ASSERT_EQ(4, row.size());

    EXPECT_FLOAT_EQ(0, matrix[0][0]);
    EXPECT_FLOAT_EQ(0.222222, matrix[0][1]);
    EXPECT_FLOAT_EQ(0.388889, matrix[3][2]);
    EXPECT_FLOAT_EQ(1, matrix[6][3]);
    EXPECT_FLOAT_EQ(1, matrix[9][1]);
}

TEST(ScoreMatrix, log_adjusted_likelihood_ratio)
{
    const unsigned number_of_sites = 18;
    detail::PWM pwm { /*number_of_sites =*/ number_of_sites };
    pwm.matrix = { {.25, .25, .25, .25},
                   {0, 0, 1, 0} };
    const double number_of_pseudo_sites=.1;

    const auto min_max = detail::log_adjusted_likelihood_ratio(pwm, uniform_bg);
    EXPECT_EQ(number_of_sites, pwm.number_of_sites);
    ASSERT_EQ(2, pwm.matrix.size());
    ASSERT_EQ(4, pwm.matrix[0].size());

    const double zero = std::log2(number_of_pseudo_sites * uniform_bg[0] / ( number_of_sites + number_of_pseudo_sites)/uniform_bg[0]);
    const double one = std::log2((number_of_sites+number_of_pseudo_sites*uniform_bg[0]) / (number_of_sites + number_of_pseudo_sites)/uniform_bg[0]);
    const double quarter = 0; // matching a base for a position where all bases are equally likely scores zero points
    EXPECT_FLOAT_EQ(quarter, pwm.matrix[0][0]);
    EXPECT_FLOAT_EQ(quarter, pwm.matrix[0][1]);
    EXPECT_FLOAT_EQ(quarter, pwm.matrix[0][2]);
    EXPECT_FLOAT_EQ(quarter, pwm.matrix[0][3]);
    EXPECT_FLOAT_EQ(zero,    pwm.matrix[1][0]);
    EXPECT_FLOAT_EQ(zero,    pwm.matrix[1][1]);
    EXPECT_FLOAT_EQ(one,     pwm.matrix[1][2]);
    EXPECT_FLOAT_EQ(zero,    pwm.matrix[1][3]);

    EXPECT_FLOAT_EQ(min_max.first,  zero);
    EXPECT_FLOAT_EQ(min_max.second, one);
}

TEST(ScoreMatrix, scale)
{
    detail::PWM pwm { /*number_of_sites =*/ 10 };
    pwm.matrix = { {0, 0, 0, 0},
                   {-8, -8, 2, -8} };
    const detail::ScaledPWM scaled = detail::scale(pwm, std::make_pair(-8.0, 2.0), 30);

    const auto& matrix = scaled.matrix;
    ASSERT_EQ(2, matrix.size());

    EXPECT_EQ(10, scaled.number_of_sites);
    EXPECT_EQ(-8, scaled.min_before_scaling);
    EXPECT_EQ(3, scaled.scale); // max - min = 10, 10*3 = 30, so scale is 3
    EXPECT_EQ(30, scaled.range);

    EXPECT_EQ(24, matrix[0][0]); // (0 - -8) * 3 = 24
    EXPECT_EQ(24, matrix[0][1]);
    EXPECT_EQ(24, matrix[0][2]);
    EXPECT_EQ(24, matrix[0][3]);

    EXPECT_EQ(0, matrix[1][0]); // (-8 - -8) * 3
    EXPECT_EQ(0, matrix[1][1]);
    EXPECT_EQ(30, matrix[1][2]); // (2 - -8) * 3
    EXPECT_EQ(0, matrix[1][3]);
}

TEST(ScoreMatrix, scaled_score)
{
    const std::vector<std::array<unsigned, AlphabetSize>> matrix =
    //  A   C   G   T
    { { 24, 24, 24, 0 },
      { 0,  0,  30, 0 } };

    EXPECT_EQ(0, detail::score(matrix, "", 0, 0));
    EXPECT_EQ(0, detail::score(matrix, "AA", 0, 0));
    EXPECT_EQ(0, detail::score(matrix, "AG", 2, 2));

    EXPECT_EQ(24, detail::score(matrix, "A", 0, 1));
    EXPECT_EQ(0,  detail::score(matrix, "T", 0, 1));
    EXPECT_EQ(0, detail::score(matrix, "N", 0, 1));
    EXPECT_EQ(0, detail::score(matrix, "Z", 0, 1));

    EXPECT_EQ(24, detail::score(matrix, "AA", 0, 2));
    EXPECT_EQ(24, detail::score(matrix, "AA", 1, 2));
    EXPECT_EQ(54, detail::score(matrix, "AG", 0, 2));
    EXPECT_EQ(54, detail::score(matrix, "ag", 0, 2));
    EXPECT_EQ(54, detail::score(matrix, "AGN", 0, 2));
    EXPECT_EQ(54, detail::score(matrix, "NAGN", 1, 3));
}

TEST(ScoreMatrix, probability_distribution)
{
    using namespace detail;

    // score of 0 is 100% probable when empty matrix
    const std::vector<std::array<unsigned, AlphabetSize>> empty_matrix;
    std::vector<double> probabilities = probability_distribution(empty_matrix, uniform_bg);
    ASSERT_EQ(1, probabilities.size());
    EXPECT_FLOAT_EQ(1, probabilities[0]);

    // score of 0 is 100% probable for a zero matrix
    const std::vector<std::array<unsigned, AlphabetSize>> zero_matrix =
    //  A   C   G   T
    { { 0,  0,  0, 0 },
      { 0,  0,  0, 0 } };
    probabilities = probability_distribution(zero_matrix, uniform_bg);
    ASSERT_EQ(1, probabilities.size());
    EXPECT_FLOAT_EQ(1, probabilities[0]);

    const std::vector<std::array<unsigned, AlphabetSize>> length_one_matrix =
    //  A   C   G   T
    { { 0,  0,  1, 0 } };

    // sequence length 1 with max 1 per base can have values 0 or 1
    // value of 0 with 75% probability
    // value of 1 with 25% probability
    probabilities = probability_distribution(length_one_matrix, uniform_bg);
    ASSERT_EQ(2, probabilities.size());
    EXPECT_FLOAT_EQ(.75, probabilities[0]);
    EXPECT_FLOAT_EQ(.25, probabilities[1]);

    const std::vector<std::array<unsigned, AlphabetSize>> length_two_matrix =
    //  A   C   G   T
    { { 0,  0,  1, 1 },
      { 1,  0,  1, 0 } };
    // Scores for every possible sequence:
    // AA: 1, AC: 0, AG: 1, AT: 0
    // CA: 1, CC: 0, CG: 1, CT: 0
    // GA: 2, GC: 1, GG: 2, GT: 1
    // TA: 2, TC: 1, TG: 2, TT: 1
    // 16 sequences
    // 4 ways to score 0: 25%
    // 8 ways to score 1: 50%
    // 4 ways to score 2: 25%
    probabilities = probability_distribution(length_two_matrix, uniform_bg);
    ASSERT_EQ(3, probabilities.size());
    EXPECT_FLOAT_EQ(.25, probabilities[0]);
    EXPECT_FLOAT_EQ(.50, probabilities[1]);
    EXPECT_FLOAT_EQ(.25, probabilities[2]);
}

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

/* The MIT License (MIT) 

   Copyright (c) 2015 John DiMatteo (jdimatteo@gmail.com)

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE. 
 */