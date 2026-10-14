#pragma once

#include <cmath>

namespace clq {

/**
 @brief  Sum of doubles with Neumaier compensation (an improved Kahan summation).

 The error of the result does not increase with the number of terms, and does
 not depend on the order of the terms. Thus the quality is more equal between
 platforms and compilers. Neumaier's version is also correct when a term is
 larger than the sum, which occurs here because gain and loss terms have
 opposite signs.

 Do not compile with -ffast-math: it can remove the compensation.
 */
class CompensatedSum {
public:
    CompensatedSum() : sum_(0.0), compensation_(0.0) {}

    void add(double x)
    {
        double t = sum_ + x;
        if (std::fabs(sum_) >= std::fabs(x)) {
            compensation_ += (sum_ - t) + x;
        } else {
            compensation_ += (x - t) + sum_;
        }
        sum_ = t;
    }

    double value() const { return sum_ + compensation_; }

private:
    double sum_;
    double compensation_;
};

}
