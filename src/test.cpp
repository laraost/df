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

// Some example trajectories at distance 1.
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

// A slightly more complex example.
bool test_2() {
    std::vector<point_t> trajectory_p = {{0,0}, {1,0}, {2,0}, {3,0}};
    std::vector<point_t> trajectory_q = {{1.2,1.1}, {1,1.1}, {2.1,1}, {3.1,1}, {4,0}};
    auto frechet = df::compute_discrete_frechet<squared_euclid>(trajectory_p.begin(), trajectory_p.end(),
                                                                trajectory_q.begin(), trajectory_q.begin() + 2);  
    const auto true_distance = squared_euclid()(trajectory_q[1], trajectory_p[3]);
    if (frechet != true_distance) {
        std::cout << "Test 2: expected Fréchet distance " << true_distance << ".\n";
        return false;
    }
    return true;
}

// Takes a lambda as distance function.
bool test_3() {
    std::vector<point_t> trajectory_p = {{0,0}, {1,0}, {2,0}, {3,0}};
    std::vector<point_t> trajectory_q = {{1.2,1.1}, {1,1.1}, {2.1,1}, {3.1,1}, {4,0}};
    auto frechet = df::compute_discrete_frechet(trajectory_p.begin(), trajectory_p.end(),
                                                trajectory_q.begin(), trajectory_q.begin() + 2,
                                                [](const point_t &, const point_t &) { return 0; } );  
    if (frechet != 0) {
        std::cout << "Test 3: expected Fréchet distance 0.\n";
        return false;
    }
    return true;
}

// One trajectory is a single point
bool test_4() {
    bool result = true;
    std::vector<point_t> trajectory_p = {{0,0}, {1,0}, {2,0}, {3,0}};
    std::vector<point_t> trajectory_q = {{1.2,1.1}, {1,1.1}, {2.1,1}, {3.1,1}, {4,0}};
    auto frechet_1 = df::compute_discrete_frechet<squared_euclid>(trajectory_p.begin(), ++trajectory_p.begin(),
                                                                trajectory_q.begin(), trajectory_q.end());  
    const auto true_distance_1 = squared_euclid()(trajectory_p[0], trajectory_q[4]);
    if (frechet_1 != true_distance_1) {
        std::cout << "Test 4: expected Fréchet distance " << true_distance_1 << ".\n";
        result = false;
    }
    auto frechet_2 = df::compute_discrete_frechet<squared_euclid>(trajectory_p.begin(), trajectory_p.end(),
                                                                  ++trajectory_q.begin(), trajectory_q.begin() + 2);  
    const auto true_distance_2 = squared_euclid()(trajectory_p[3], trajectory_q[1]);
    if (frechet_2 != true_distance_2) {
        std::cout << "Test 4: expected Fréchet distance " << true_distance_2 << ".\n";
        result = false;
    }
    return result;
}

int main() {
    int errors = 0;
    errors += static_cast<int>(!test_1());
    errors += static_cast<int>(!test_2());
    errors += static_cast<int>(!test_3());
    errors += static_cast<int>(!test_4());

    std::cout << "Errors: " << errors << "\n";

    return errors != 0;
}
