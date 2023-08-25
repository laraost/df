#include "discrete_frechet.h"

#include <cassert>
#include <iostream>
#include <utility>
#include <vector>

using point_t = std::pair<double, double>;

struct squared_euclid {

    double operator()(const point_t &a, const point_t &b) const {
        return (a.first - b.first)*(a.first - b.first) + (a.second - b.second)*(a.second - b.second);
    }

};

bool test_1() {
    std::vector<point_t> trajectory_p = {{0,0}, {1,0}, {2,0}, {3,0}};
    std::vector<point_t> trajectory_q = {{0,1}, {1,1}, {2,1}, {3,1}, {4,0}};
    auto frechet = df::compute_discrete_frechet<squared_euclid>(trajectory_p.begin(), trajectory_p.end(),
                                                                trajectory_q.begin(), trajectory_q.end());  
    if (frechet != 1) {
        std::cout << "Test 1: expected Fréchet distance 1.\n";
        return false;
    }
    return true;
}

bool test_2() {
    std::vector<point_t> trajectory_p = {{0,0}, {1,0}, {2,0}, {3,0}};
    std::vector<point_t> trajectory_q = {{1.2,1.1}, {1,1.1}, {2.1,1}, {3.1,1}, {4,0}};
    auto frechet = df::compute_discrete_frechet<squared_euclid>(trajectory_p.begin(), trajectory_p.end(),
                                                                trajectory_q.begin(), trajectory_q.begin() + 2);  
    const auto true_distance = squared_euclid()(trajectory_q[1], trajectory_p[3]);
    if (frechet != true_distance) {
        std::cout << "Test 1: expected Fréchet distance " << true_distance << ".\n";
        return false;
    }
    return true;
}

int main() {
    int errors = 0;
    errors += static_cast<int>(!test_1());
    errors += static_cast<int>(!test_2());

    std::cout << "Errors: " << errors << "\n";

    return errors != 0;
}
