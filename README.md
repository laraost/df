# Discrete Fréchet Distance

This is a header-only implementation of the dynamic programming algorithm by Eiter and Mannila (1994) to compute the discrete Fréchet distance.
Requires at least C++14.

## Usage

Place the header `src/discrete_frechet.h` at a convenient place.
The function to compute the Fréchet distance takes as input two trajectories, defined by a pair of (forward) iterators each,
plus a distance function.

Example (see also `src/test.cpp`):
```
using point_t = std::pair<double, double>;
std::vector<point_t> trajectory_p = {{0,0}, {1,0}, {2,0}, {3,0}};
std::vector<point_t> trajectory_q = {{0,1}, {1,1}, {2,1}, {3,1}, {4,0}};
auto frechet = df::compute_discrete_frechet<squared_euclid>(trajectory_p.begin(), trajectory_p.end(),
                                                            trajectory_q.begin(), trajectory_q.end(),
                                                            [](const point_t &a, const point_t &b) {
                                                                return (a.first - b.first)*(a.first - b.first) +
                                                                    (a.second - b.second)*(a.second - b.second);
                                                            });  
```
