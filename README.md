# Unscented Kalman Filter Project
The Unscented Kalman Filter (UKF) is an improvement over the Extended Kalman Filter (EKF) in the way it handles non-linear motion. The EKF predicts non-linear motion by calculating a linear approximation of a slice of it (the Jacobian), which might not work well on highly non-linear motion models. The UKF tries to solve this problem by sampling several points (sigma points) around the mean of the measurements and feeding them through the non-linear motion model, and then reconstructing the new mean and covariance based on these transformed sigma points.

---

## Dependencies

* cmake >= v3.5
* make >= v4.1
* gcc/g++ >= v5.4

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF path/to/input.txt path/to/output.txt`. You can find
   some sample inputs in 'data/'.
    - eg. `./UnscentedKF ../data/obj_pose-laser-radar-synthetic-input.txt`
