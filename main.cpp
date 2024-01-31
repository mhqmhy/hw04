#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>

#define N 48
float frand() {
    return (float)rand() / RAND_MAX * 2 - 1;
}

struct Star {
    float px[N], py[N], pz[N];
    float vx[N], vy[N], vz[N];
    float mass[N];
    float padding[512 - 336];
    // 48 * 7 = 336
};

// std::vector<Star> stars;
Star stars;

void init() {
    for (std::size_t i = 0; i < 48; i++) {
        stars.px[i] = frand();
        stars.py[i] = frand();
        stars.pz[i] = frand();
        stars.vx[i] = frand();
        stars.vy[i] = frand();
        stars.vz[i] = frand();
        stars.mass[i] = frand() + 1;
        // stars.push_back({
        //     frand(), frand(), frand(),
        //     frand(), frand(), frand(),
        //     frand() + 1,
        // });
    }
}

 float G = 0.001;
 float eps = 0.001;
 float dt = 0.01;
 float eps2 = eps * eps;
 float G_dt = G * dt;

void step() {
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC unroll 8
#elif defined(_MSC_VER)
#pragma unroll 8
#endif
    // for (auto &star: stars) {
    //     for (auto &other: stars) {
    //         float dx = other.px - star.px;
    //         float dy = other.py - star.py;
    //         float dz = other.pz - star.pz;
    //         float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
    //         d2 *= sqrt(d2);
    //         star.vx += dx * other.mass * G * dt / d2;
    //         star.vy += dy * other.mass * G * dt / d2;
    //         star.vz += dz * other.mass * G * dt / d2;
    //     }
    // }
    for (std::size_t i = 0; i < 48; i++) {
        for (std::size_t j = 0; j < 48; j++) {
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps2;
            d2 *= std::sqrt(d2);
            stars.vx[i] += dx * stars.mass[j] * G_dt / d2;
            stars.vy[i] += dy * stars.mass[j] * G_dt / d2;
            stars.vz[i] += dz * stars.mass[j] * G_dt / d2;
        }
    }

#if defined(__GNUC__) || defined(__clang__)
#pragma GCC unroll 8
#elif defined(_MSC_VER)
#pragma unroll 8
#endif
    for (std::size_t i = 0; i < N; i++) {
        stars.px[i] += stars.vx[i] * dt;
        stars.py[i] += stars.vy[i] * dt;
        stars.pz[i] += stars.vz[i] * dt;
    }
    // for (auto &star: stars) {
    //     star.px += star.vx * dt;
    //     star.py += star.vy * dt;
    //     star.pz += star.vz * dt;
    // }
}

float calc() {
    float energy = 0;
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC unroll 8
#elif defined(_MSC_VER)
#pragma unroll 8
#endif
    // for (auto &star: stars) {
    //     float v2 = star.vx * star.vx + star.vy * star.vy + star.vz * star.vz;
    //     energy += star.mass * v2 / 2;
    //     for (auto &other: stars) {
    //         float dx = other.px - star.px;
    //         float dy = other.py - star.py;
    //         float dz = other.pz - star.pz;
    //         float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
    //         energy -= other.mass * star.mass * G / sqrt(d2) / 2;
    //     }
    // }

    for (std::size_t i = 0; i < N; i++) {
        float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
        energy += stars.mass[i] * v2 / 2;
        for (std::size_t j = 0; j < N; j++) {
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps2;
            energy -= stars.mass[j] * stars.mass[i] * G / std::sqrt(d2) / 2;

        }
    }
    return energy;
}

template <class Func>
long benchmark(Func const &func) {
    auto t0 = std::chrono::steady_clock::now();
    func();
    auto t1 = std::chrono::steady_clock::now();
    auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
    return dt.count();
}

int main() {
    init();
    printf("Initial energy: %f\n", calc());
    auto dt = benchmark([&] {
        for (int i = 0; i < 100000; i++)
            step();
    });
    printf("Final energy: %f\n", calc());
    printf("Time elapsed: %ld ms\n", dt);
    return 0;
}

/*
Initial energy: -13.414000
Final energy: -13.356842
Time elapsed: 2015 ms
*/