//
// Created by phili on 23.08.2024.
//

#include <chrono>
#include <iostream>
#include <random>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <filesystem>
#include <valarray>

// Step 1: Insert custom values here. There are no constraints.

// Dimensions [mm]
double a = 800.0;
double b = 400.0;
double plyThickness = 0.184;

// Force per width [kN/mm]
double Nx = -0.2 * 10 * 10 * 10;
double Ny = -0.01 * 10 * 10 * 10;
double Tau = 0.1 * 10 * 10 * 10;

// Material Properties [MPa]
double E1 = 130000;
double E2 = 10000;
double G12 = 5000;

double nu12 = 0.3;

// Constants
double PI = 3.14159;

// Optimization Parameters
double RF_weight = 1;
double size_weight = 100000;
int initial_population = 15;
int initial_max_stack_size = 50;
bool fillPopulation = true;

// Necessary Structs
struct QMatrix {
    double Q11;
    double Q22;
    double Q12;
    double Q66;
};

struct QBarMatrix {
    double Q11 = 0;
    double Q22 = 0;
    double Q66 = 0;
    double Q12 = 0;
    double Q16 = 0;
    double Q26 = 0;

    QBarMatrix() = default;

    QBarMatrix(const QMatrix &qMatrix, const double angle) {
        const double rad = angle / 180 * PI;
        const double c = std::cos(rad);
        const double s = std::sin(rad);
        const double c2 = c * c;
        const double s2 = s * s;
        const double c4 = c2 * c2;
        const double s4 = s2 * s2;
        const double c2s2 = c2 * s2;

        Q11 = qMatrix.Q11 * c4 + 2 * (qMatrix.Q12 + 2 * qMatrix.Q66) * c2s2 + qMatrix.Q22 * s4;
        Q22 = qMatrix.Q11 * s4 + 2 * (qMatrix.Q12 + 2 * qMatrix.Q66) * c2s2 + qMatrix.Q22 * c4;
        Q12 = (qMatrix.Q11 + qMatrix.Q22 - 4 * qMatrix.Q66) * c2s2 + qMatrix.Q12 * (c4 + s4);
        Q66 = (qMatrix.Q11 + qMatrix.Q22 - 2 * qMatrix.Q12 - 2 * qMatrix.Q66) * c2s2 + qMatrix.Q66 * (c4 + s4);
        Q16 = (qMatrix.Q11 - qMatrix.Q12 - 2 * qMatrix.Q66) * c * c2 * s + (qMatrix.Q12 - qMatrix.Q22 + 2 * qMatrix.Q66)
              * s * s2 * c;
        Q26 = (qMatrix.Q11 - qMatrix.Q12 - 2 * qMatrix.Q66) * s * s2 * c + (qMatrix.Q12 - qMatrix.Q22 + 2 * qMatrix.Q66)
              * c * c2 * s;
    }
};

enum PlyOrientation {
    ANGLE_0,
    ANGLE_PLUS_45,
    ANGLE_MINUS_45,
    ANGLE_90,
    NAN_ANGLE
};

// Initial Calculations
double nu21;

QMatrix Q;

std::unordered_map<PlyOrientation, QBarMatrix> QBarMatrices;

std::unordered_map<PlyOrientation, int> angles;

struct DMatrix {
    double D11 = 0;
    double D22 = 0;
    double D12 = 0;
    double D66 = 0;
    double D16 = 0;
    double D26 = 0;

    void updateDMatrix(double z1, double z2, QBarMatrix Q) {
        z1 = std::pow(z1, 3);
        z2 = std::pow(z2, 3);
        D11 += 1.0 / 3 * Q.Q11 * (z1 - z2);
        D22 += 1.0 / 3 * Q.Q22 * (z1 - z2);
        D12 += 1.0 / 3 * Q.Q12 * (z1 - z2);
        D66 += 1.0 / 3 * Q.Q66 * (z1 - z2);
        D16 += 1.0 / 3 * Q.Q16 * (z1 - z2);
        D26 += 1.0 / 3 * Q.Q26 * (z1 - z2);
    }
};

struct Stack {
    // The stacking sequence only includes the upper half of the stack due to symmetry
    std::vector<PlyOrientation> stackingSequence;
    DMatrix D;
    int stackSize = static_cast<int>(stackingSequence.size() * 2);
    double fitness{};
    double fitnessReserve = 0;
    double fitnessStackSize = 0;
    double RF = 0;
    std::vector<Stack> solutions;
    double rank = 0;
    int domination = 0;
    double crowdingDistance = 0;

public:
    double calculateFitness() {
        stackSize = stackingSequence.size();
        if (D.D11 == 0) {
            setupDMatrix();
        }
        const double RF_biax = std::abs(
            calculateSigmaCrit() / (1.5 * Nx / (2 * static_cast<double>(stackingSequence.size()) * plyThickness)));
        const double RF_shear = std::abs(
            calculateTauCrit() / (1.5 * (Tau / (2 * static_cast<double>(stackingSequence.size()) * plyThickness))));
        RF = 1 / (1 / RF_biax + 1 / RF_shear * 1 / RF_shear);
        if (std::isinf(RF)) {
            std::cout << "WARNING: Encountered 0 in denominator." << std::endl;
            fitness = -1;
            return fitness;
        }
        fitnessReserve = 1.0 / (std::abs((1 - RF)) * RF_weight);
        fitnessStackSize = 1.0 / (stackSize * size_weight);
        if (RF < 1) {
            rank = 2;
            fitness = -1;
            return fitness;
        }
        fitness = fitnessReserve * fitnessReserve;
        return fitness;
    }

private:
    void setupDMatrix() {
        // Calculate the D-Matrix
        for (int i = static_cast<int>(stackingSequence.size()); i > 0; i--) {
            double z1 = i * plyThickness;
            double z2 = (i - 1) * plyThickness;
            D.updateDMatrix(z1, z2, QBarMatrices[stackingSequence.at(stackingSequence.size() - i)]);
            z1 = -(i - 1) * plyThickness;
            z2 = -i * plyThickness;
            D.updateDMatrix(z1, z2, QBarMatrices[stackingSequence.at(stackingSequence.size() - i)]);
        }
    }

    [[nodiscard]] double calculateTauCrit() const {
        const double stiffness_ratio = std::sqrt(D.D11 * D.D22) / (D.D12 + 2 * D.D66);
        if (stiffness_ratio >= 1) {
            return 4 / (stackSize * plyThickness * b * b) * (
                       std::sqrt(D.D11 * D.D22 * D.D22 * D.D22) * std::sqrt(D.D11 * D.D22 * D.D22 * D.D22) * (
                           8.12 + 5.05 / stiffness_ratio));
        }
        return 4 / (stackSize * plyThickness * b * b) * (
                   sqrt(D.D22 * (D.D12 + 2 * D.D66)) * (
                       11.7 + 0.532 * stiffness_ratio + 0.938 * stiffness_ratio * stiffness_ratio));
    }


    [[nodiscard]] double calculateSigmaCrit() const {
        double sCritM = 0;
        double sCritN = 0;
        bool flag = false;
        int m = 1;
        int n = 1;
        double _ = sigmaCrit(2, 1);
        while (true) {
            double _ = sigmaCrit(m, 1);
            if (!flag) {
                sCritM = _;
                flag = true;
            } else {
                if (_ < std::abs(sCritM)) {
                    sCritM = _;
                } else break;
            }
            m++;
        }
        flag = false;
        while (true) {
            double _ = sigmaCrit(1, n);
            if (!flag) {
                sCritN = _;
                flag = true;
            } else {
                if (std::abs(_) < std::abs(sCritN)) {
                    sCritN = _;
                } else break;
            }
            n++;
        }
        if (sCritM < sCritN) {
            return sCritM;
        }
        return sCritN;
    }

    [[nodiscard]] double sigmaCrit(const int m, const int n) const {
        double alpha = a / b;
        double beta = Ny / Nx;
        return PI * PI / (b * b * stackSize * plyThickness) * 1 / (
                   (m / alpha) * (m / alpha) + beta * n * n) * (
                   D.D11 * (m / alpha) * (m / alpha) * (m / alpha) * (m / alpha) + 2 * (D.D12 + D.D66) * (m * n / alpha)
                   * (m * n / alpha) + D.D22 * n * n * n * n);
    }
};

void initialCalculations() {
    // Initial Calculations
    nu21 = nu12 * E2 / E1;

    Q = {E1 / (1 - nu12 * nu21), E2 / (1 - nu12 * nu21), nu12 * E2 / (1 - nu12 * nu21), G12};

    QBarMatrices = {
        {ANGLE_0, QBarMatrix(Q, 0)},
        {ANGLE_PLUS_45, QBarMatrix(Q, 45)},
        {ANGLE_MINUS_45, QBarMatrix(Q, -45)},
        {ANGLE_90, QBarMatrix(Q, 90)}
    };

    angles = {
        {ANGLE_0, 0},
        {ANGLE_PLUS_45, 45},
        {ANGLE_MINUS_45, -45},
        {ANGLE_90, 90}
    };
}

bool checkStackValidity(const std::vector<PlyOrientation> &stack) {
    std::unordered_map<PlyOrientation, double> orientationSize;
    int stackSize = static_cast<double>(stack.size());
    for (PlyOrientation ply: stack) {
        orientationSize[ply] += 1.0 / (2 * stackSize);
    }
    for (const auto &size: orientationSize) {
        if (size.second < 0.1) {
            // std::cout << "Stack has a ply with ply share strictly less than 10%. It will be ignored." << std::endl;
            return false;
        }
    }
    return true;
}

Stack generateStack(int size) {
    /*
     * Input: number of plies in the upper half of the stack
     * Return: stack following given boundary conditions
     * Boundary Conditions:
     *  1. symmetric
     *  2. balanced
     *  3. min. 10% each ply share
     *  4. given ply thickness
     */
    // Check edge cases
    if (size <= 0) throw std::runtime_error("Size of stack must be strictly larger than 0");

    std::vector<PlyOrientation> stack;
    const std::vector<PlyOrientation> orientations = {ANGLE_0, ANGLE_PLUS_45, ANGLE_MINUS_45, ANGLE_90};

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, static_cast<int>(orientations.size() - 1));

    for (int i = 0; i < size; i++) {
        stack.push_back(orientations[dist(gen)]);
    }
    if (checkStackValidity(stack)) {
        Stack s = {stack};
        s.calculateFitness();
        return s;
    }
    Stack s = {std::vector<PlyOrientation>(size, orientations[dist(gen)])};
    s.calculateFitness();
    return s;
}

std::vector<Stack> generatePopulation(const int size, const int maxStackSize, const int minStackSize = 1) {
    std::vector<Stack> population;

    std::random_device dev;
    std::mt19937 gen(dev());
    std::uniform_int_distribution<> dist(minStackSize, maxStackSize);

    for (int i = 0; i < size; i++) {
        population.push_back(generateStack(dist(gen)));
    }

    return population;
}

std::pair<Stack, Stack> crossover(const Stack &father, const Stack &mother) {
    const auto stackSizeFather = static_cast<size_t>(father.stackingSequence.size());
    const auto stackSizeMother = static_cast<size_t>(mother.stackingSequence.size());
    auto stackSize = 0;
    if (stackSizeFather < stackSizeMother) {
        stackSize = stackSizeFather;
    } else
        stackSize = stackSizeMother;
    Stack daughter, son;

    auto crossover_point = std::rand() % stackSize;

    for (size_t i = 0; i < stackSize; i++) {
        if (i <= crossover_point) {
            son.stackingSequence.push_back(father.stackingSequence[i]);
            daughter.stackingSequence.push_back(mother.stackingSequence[i]);
        } else {
            son.stackingSequence.push_back(mother.stackingSequence[i]);
            daughter.stackingSequence.push_back(father.stackingSequence[i]);
        }
    }

    son.calculateFitness();
    daughter.calculateFitness();

    return {son, daughter};
}

std::string plyOrientationToString(PlyOrientation ply_orientation) {
    switch (ply_orientation) {
        case ANGLE_0: return "0°";
        case ANGLE_PLUS_45: return "45°";
        case ANGLE_MINUS_45: return "-45°";
        case ANGLE_90: return "90°";
    }
    return "";
}

bool isDominant(const Stack &lhs, const Stack &rhs) {
    // Checks if lhs is at least as good as rhs or at least in one objective strictly better
    return lhs.fitnessStackSize >= rhs.fitnessStackSize && lhs.fitnessReserve >= rhs.fitnessReserve && (
               lhs.fitnessStackSize > rhs.fitnessStackSize || lhs.fitnessReserve > rhs.fitnessReserve);
}

std::vector<std::vector<Stack> > fast_non_dominated_sort(std::vector<Stack> &population) {
    std::vector<std::vector<Stack> > Fronts = {{}};
    for (int p = 0; p < population.size(); p++) {
        population[p].solutions = {};
        for (int q = 0; q < population.size(); q++) {
            if (p == q) continue;
            if (isDominant(population[p], population[q])) {
                population[p].solutions.push_back(population[q]);
            } else if (isDominant(population[q], population[p])) {
                population[p].domination++;
            }
        }
        if (population[p].domination == 0) {
            // p belongs to the first front
            population[p].rank = 0;
            population[p].calculateFitness();
            Fronts[0].push_back(population[p]);
        }
    }
    int front = 0;
    while (!Fronts[front].empty()) {
        std::vector<Stack> Q = {};
        for (int p = 0; p < Fronts[front].size(); p++) {
            for (Stack q: population[p].solutions) {
                q.domination--;
                if (q.domination == 0) {
                    q.rank = front + 1;
                    Q.push_back(q);
                }
            }
        }
        front++;
        Fronts.push_back(Q);
    }
    return Fronts;
}

void crowding_distance_assignment(std::vector<Stack> population) {
    const int l = population.size();

    std::sort(population.begin(), population.end(),
              [](const Stack &lhs, const Stack &rhs) {
                  return lhs.fitnessStackSize > rhs.fitnessStackSize;
              });
    population[0].crowdingDistance = population[l - 1].crowdingDistance = std::numeric_limits<double>::infinity();
    for (int i = 1; i < l - 1; i++) {
        population[i].crowdingDistance += (population[i + 1].fitnessStackSize - population[i - 1].fitnessStackSize) / (
            population[0].fitnessStackSize - population[l - 1].fitnessStackSize);
    }

    std::sort(population.begin(), population.end(),
              [](const Stack &lhs, const Stack &rhs) {
                  return lhs.fitnessReserve > rhs.fitnessReserve;
              });
    population[0].crowdingDistance = population[l - 1].crowdingDistance = std::numeric_limits<double>::infinity();
    for (int i = 1; i < l - 1; i++) {
        population[i].crowdingDistance += (population[i + 1].fitnessReserve - population[i - 1].fitnessReserve) / (
            population[0].fitnessReserve - population[l - 1].fitnessReserve);
    }
}

bool crowded_comparison_operator(const Stack &i, const Stack &j) {
    return (i.rank < j.rank) || ((i.rank == j.rank && i.crowdingDistance > j.crowdingDistance));
}

Stack tournament_selection(std::vector<Stack> population) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(1, static_cast<int>(population.size() - 1));

    Stack contastent1 = population[dist(gen)];
    Stack contastent2 = population[dist(gen)];

    if (crowded_comparison_operator(contastent1, contastent2)) return contastent1;
    return contastent2;
}

PlyOrientation getAngle(int angle) {
    for (const auto &pair: angles) {
        if (pair.second == angle) {
            return pair.first;
        }
    }
    return NAN_ANGLE;
}

void mutation(Stack &offspring, const int mutation = 1) {
    std::random_device rd;
    std::uniform_real_distribution<> dist(0.1, 1);
    std::uniform_int_distribution<> angle(-mutation, mutation);
    int size = offspring.stackingSequence.size() * dist(rd);
    for (int j = dist(rd); j < size; j++) {
        int newAngle = angles.at(offspring.stackingSequence[j]) + angle(rd) * 45;
        if (std::abs(newAngle) > 90) {
            newAngle = 0;
        }
        if (newAngle == -90) {
            newAngle = 90;
        }
        offspring.stackingSequence[j] = getAngle(newAngle);
    }
    size -= dist(rd) * size;
    offspring.calculateFitness();
    if (offspring.RF > 1) {
        offspring.stackingSequence.erase(offspring.stackingSequence.end() - size, offspring.stackingSequence.end());
        return;
    }
    Stack addon = generateStack(size + 1);
    offspring.stackingSequence.insert(offspring.stackingSequence.end(), addon.stackingSequence.begin(),
                                      addon.stackingSequence.end());
}

std::vector<Stack> make_new_pop(std::vector<Stack> population) {
    std::vector<Stack> new_population;
    int N = population.size();

    while (new_population.size() < N) {
        Stack father = tournament_selection(population);
        Stack mother = tournament_selection(population);

        auto [son, daughter] = crossover(father, mother);

        mutation(son);
        mutation(daughter);

        new_population.push_back(son);
        if (new_population.size() < N) new_population.push_back(daughter);
    }
    return new_population;
}

void NSGA2(std::vector<Stack> &Pt, int N, int generations) {
    std::vector<Stack> Qt;
    int t = 0;
    while (t < generations) {
        if (t % 1 == 0)
            std::cout << t << std::endl;
        std::vector<Stack> Rt = Pt;
        Rt.insert(Rt.end(), Qt.begin(), Qt.end());

        std::vector<std::vector<Stack>> F = fast_non_dominated_sort(Rt); // non-dominant fronts

        std::vector<Stack> Pt_1 = {};
        int i = 0;
        while (i < F.size() && Pt_1.size() + F[i].size() <= N) {
            if (F[i].empty()) {
                i--;
                break;
            }
            crowding_distance_assignment(F[i]);
            Pt_1.insert(Pt_1.end(), F[i].begin(), F[i].end());
            i++;
        }
        std::sort(F[i].begin(), F[i].end(), crowded_comparison_operator);
        if (F[i].size() >= N - Pt_1.size())
            Pt_1.insert(Pt_1.end(), F[i].begin(), F[i].begin() + (N - Pt_1.size()));

        Qt = make_new_pop(Pt_1);

        Pt = Pt_1;

        t++;
    }
}


int main() {
    const auto start = std::chrono::high_resolution_clock::now();
    initialCalculations();

    std::cout << "Step 1: Generating population..." << std::endl;
    std::vector<Stack> population = generatePopulation(initial_population, initial_max_stack_size);

    NSGA2(population, population.size(), 5000);

    // Finished
    auto stop = std::chrono::high_resolution_clock::now();
    // if (converged) {
    //     std::cout << "--- Converged! ---" << std::endl;
    // } else {
    //     std::cout << "--- Optimizer finished before completion. ---" << std::endl;
    // }

    std::cout << "Number of solutions: " << population.size() << std::endl;
    std::cout << "Found following Solution: " << std::endl;

    std::sort(population.begin(), population.end(),
              [](const Stack &lhs, const Stack &rhs) {
                  return lhs.rank < rhs.rank;
              });

    population[0].calculateFitness();
    std::cout << "  RF: " << population[0].RF << std::endl;
    std::cout << "  Rank: " << population[0].rank << std::endl;
    std::cout << "  Stack Size: " << 2 * population[0].stackingSequence.size() << std::endl;
    std::cout << "  Stacking Sequence: " << std::endl;
    for (const PlyOrientation ply: population[0].stackingSequence) {
        std::string plyOrientation = plyOrientationToString(ply);
        std::cout << "      " << plyOrientation << std::endl;
    }
    std::cout << "  --- Mid Plane --- " << std::endl;
    // for (Stack s: population) {
    //     s.calculateFitness();
    //     std::cout << "  RF: " << s.RF << std::endl;
    //     std::cout << "  Fitness: " << s.rank << std::endl;
    //     std::cout << "  Stack Size: " << 2 * s.stackingSequence.size() << std::endl;
    //     std::cout << "  Stacking Sequence: " << std::endl;
    //     for (const PlyOrientation ply: s.stackingSequence) {
    //         std::string plyOrientation = plyOrientationToString(ply);
    //         std::cout << "      " << plyOrientation << std::endl;
    //     }
    //     std::cout << "  --- Mid Plane --- " << std::endl;
    // }

    std::cout << "\n\nTime taken: " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() <<
            "ms" << std::endl;

    return 0;
}
