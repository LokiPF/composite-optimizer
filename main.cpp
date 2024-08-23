//
// Created by phili on 23.08.2024.
//

#include <iostream>
#include <random>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <valarray>

double a = 400.0;
double b = 150.0;

double Nx = -1200;
double Ny = -800;
double Tau = 800;

double plyThickness = 0.184;

double E1 = 130000;
double E2 = 10000;
double G12 = 5000;

double nu12 = 0.3;

double PI = 3.14159;

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
    ANGLE_90
};

// Initial Calculations
double nu21;

QMatrix Q;

std::pmr::unordered_map<PlyOrientation, QBarMatrix> QBarMatrices;

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
    std::pmr::vector<PlyOrientation> stackingSequence;
    DMatrix D;
    double fitness{};
    double fitnessReserve = 0;
    double fitnessStackSize = 0;

public:
    double calculateFitness() {
        if (D.D11 == 0) {
            setupDMatrix();
        }
        const double RF_biax = std::abs(
            calculateSigmaCrit() / (1.5 * Nx / (static_cast<double>(stackingSequence.size()) * plyThickness)));
        const double RF_shear = std::abs(calculateTauCrit() / (1.5 * Tau));
        const double RF = 1 / (1 / RF_biax + 1 / RF_shear * 1 / RF_shear);
        fitnessReserve = 1.0 / std::abs((1 - RF));
        fitnessStackSize = 1.0 / (stackingSequence.size() * 100000);
        if (RF < 1) {
            fitnessReserve = 0;
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

    double calculateTauCrit() const {
        const double stiffness_ratio = std::sqrt(D.D11 * D.D22) / (D.D12 + 2 * D.D66);
        if (stiffness_ratio >= 1) {
            return 4 / (static_cast<double>(stackingSequence.size()) * plyThickness * b * b) * (
                       std::sqrt(D.D11 * D.D22 * D.D22 * D.D22) * std::sqrt(D.D11 * D.D22 * D.D22 * D.D22) * (
                           8.12 + 5.05 / stiffness_ratio));
        }
        return 4 / (static_cast<double>(stackingSequence.size()) * plyThickness * b * b) * (
                   sqrt(D.D22 * (D.D12 + 2 * D.D66)) * (
                       11.7 + 0.532 * stiffness_ratio + 0.938 * stiffness_ratio * stiffness_ratio));
    }


    double calculateSigmaCrit() {
        double sCrit = 0;
        bool flag = false;
        int m = 1;
        int n = 1;
        while (true) {
            double _ = sigmaCrit(m, 1);
            if (!flag) {
                sCrit = _;
                flag = true;
            } else {
                if (_ < std::abs(sCrit)) {
                    sCrit = _;
                } else break;
            }
            m++;
        }
        flag = false;
        while (true) {
            double _ = sigmaCrit(1, n);
            if (!flag) {
                sCrit = _;
                flag = true;
            } else {
                if (std::abs(_) < std::abs(sCrit)) {
                    sCrit = _;
                } else break;
            }
            n++;
        }
        return sCrit;
    }

    [[nodiscard]] double sigmaCrit(const int m, const int n) const {
        double alpha = a / b;
        double beta = Ny / Nx;
        return PI * PI / (b * b * static_cast<double>(stackingSequence.size()) * plyThickness) * 1 / (
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
}

bool checkStackValidity(const std::pmr::vector<PlyOrientation> &stack) {
    std::pmr::unordered_map<PlyOrientation, double> orientationSize;
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

    std::pmr::vector<PlyOrientation> stack;
    const std::pmr::vector<PlyOrientation> orientations = {ANGLE_0, ANGLE_PLUS_45, ANGLE_MINUS_45, ANGLE_90};

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
    Stack s = {std::pmr::vector<PlyOrientation>(size, orientations[dist(gen)])};
    s.calculateFitness();
    return s;
}

std::pmr::vector<Stack> generatePopulation(const int size, const int maxPlySize) {
    std::pmr::vector<Stack> population;

    std::random_device dev;
    std::mt19937 gen(dev());
    std::uniform_int_distribution<> dist(1, maxPlySize);

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

int main() {
    initialCalculations();

    std::cout << "OPTIMIZER STARTED..." << std::endl;

    std::cout << "Step 1: Generating population..." << std::endl;
    std::pmr::vector<Stack> population = generatePopulation(10000, 200);

    bool converged = false;

    while (!converged) {
        std::cout << "Step 2: Selecting best 15%..." << std::endl;
        std::sort(population.begin(), population.end(),
                  [](const Stack &lhs, const Stack &rhs) {
                      return lhs.fitness > rhs.fitness;
                  });
        auto top = static_cast<size_t>(0.3 * static_cast<int>(population.size()));
        std::pmr::vector<Stack> elite;
        for (int i = 0; i < top; i++) {
            elite.push_back(population[i]);
        }

        population = elite;
        if (population.size() <= 7) {
            converged = true;
            break;
        }
        std::cout << "Step 3: Starting crossover..." << std::endl;
        for (int i = 0; i < elite.size() - 1; i += 2) {
            auto [son, daughter] = crossover(elite[i], elite[i + 1]);
            population.push_back(son);
            population.push_back(daughter);
        }
    }

    std::cout << "Step 4: Converged!" << std::endl;

    std::cout << "Found following Solution: " << std::endl;

    std::cout << "  RF: " << 1 / population[0].fitnessReserve + 1 << std::endl;
    std::cout << "  Fitness: " << population[0].fitness << std::endl;
    std::cout << "  Stack Size: " << 2 * population[0].stackingSequence.size() << std::endl;
    std::cout << "  Stacking Sequence: " << std::endl;
    for (const PlyOrientation ply: population[0].stackingSequence) {
        std::string plyOrientation = plyOrientationToString(ply);
        std::cout << "      " << plyOrientation << std::endl;
    }


    return 0;
}
