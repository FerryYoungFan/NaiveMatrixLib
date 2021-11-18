#include "matBasic_real.hpp"    // Real matrices
#include "matBasic_complex.hpp" // Complex matrices

void realMatrixExamples()
{
    std::cout << "\n\n******************** Real Matrix Examples ********************\n\n";
    i_real_matrix matA = {
        {1.0, 1.0, 4.0},
        {4.0, 1.0, 5.0},
        {8.0, 1.0, 0.0}};
    showMatrix(matA, "A");
    std::cout << "rank(A) = " << rank(matA) << "\n";
    std::cout << "det(A) = " << det(matA) << "\n";
    showMatrix(inv(matA), "inv(A)");
    showMatrix(pinv(matA), "pinv(A)");
    showMatrix(pinv2(matA), "pinv2(A)");
    std::cout << "\n\n";

    i_real_matrix matb = {
        {1.0},
        {2.0},
        {3.0}};
    showMatrix(matb, "b");
    showMatrix(leftDiv(matA, matb), "A \\ b");
    showMatrix(matMul(inv(matA), matb), "inv(A) * b");
    showMatrix(matMul(pinv(matA), matb), "pinv(A) * b");
    showMatrix(matMul(pinv2(matA), matb), "pinv2(A) * b (== A \\ b)");
    std::cout << "\n\n";

    i_real_matrix matB = {
        {1.0, 1.0, 4.0},
        {4.0, 1.0, 5.0},
        {0.0, 0.0, 0.0}};
    showMatrix(matB, "B");
    std::cout << "rank(B) = " << rank(matB) << "\n";
    std::cout << "det(B) = " << det(matB) << "\n";
    // showMatrix(inv(matB), "inv(B)"); // Error: singular matrix
    // showMatrix(pinv(matB), "pinv(B)"); // Error: not a full-rank matrix
    showMatrix(pinv2(matB), "pinv2(B)");
    showMatrix(leftDiv(matB, matb), "B \\ b"); // Warning: singular matrix
    showMatrix(matMul(pinv2(matB), matb), "pinv2(B) * b (== B \\ b)"); // Warning: singular matrix
    std::cout << "\n\n";

    i_real_matrix matC = {
        {1.0, 1.0},
        {4.0, 5.0},
        {1.0, 4.0}};
    showMatrix(matC, "C");
    std::cout << "rank(C) = " << rank(matC) << "\n";
    // std::cout << "det(C) = " << det(matC) << "\n"; // Error: not a square matrix
    // showMatrix(inv(matC), "inv(C)"); // Error: not a square matrix
    showMatrix(pinv(matC), "pinv(C)");
    showMatrix(pinv2(matC), "pinv2(C)");
    showMatrix(leftDiv(matC, matb), "C \\ b");
    showMatrix(matMul(pinv(matC), matb), "pinv(C) * b");
    showMatrix(matMul(pinv2(matC), matb), "pinv2(C) * b (== C \\ b)");
    std::cout << "\n\n";
}

void complexMatrixExamples()
{
    std::cout << "\n\n******************** Complex Matrix Examples ********************\n\n";
    i_complex_matrix matA = {
        {{1.0, 0.0}, {0.0, 1.0}, {4.0, 0.0}},
        {{0.0, 5.0}, {1.0, 0.0}, {4.0, 0.0}},
        {{8.0, 0.0}, {1.0, 1.0}, {0.0, 0.0}}};
    showMatrix(matA, "A");
    std::cout << "rank(A) = " << rank(matA) << "\n";
    std::cout << "det(A) = " << det(matA) << "\n";
    showMatrix(inv(matA), "inv(A)");
    showMatrix(pinv(matA), "pinv(A)");
    showMatrix(pinv2(matA), "pinv2(A)");
    std::cout << "\n\n";

    i_complex_matrix matb = {
        {{1.0, 0.0}},
        {{1.0, 2.0}},
        {{0.0, 3.0}}};
    showMatrix(matb, "b");
    showMatrix(leftDiv(matA, matb), "A \\ b");
    showMatrix(matMul(inv(matA), matb), "inv(A) * b");
    showMatrix(matMul(pinv(matA), matb), "pinv(A) * b");
    showMatrix(matMul(pinv2(matA), matb), "pinv2(A) * b (== A \\ b)");
    std::cout << "\n\n";

    i_complex_matrix matB = {
        {{1.0, 0.0}, {0.0, 1.0}, {4.0, 0.0}},
        {{0.0, 5.0}, {1.0, 0.0}, {4.0, 0.0}},
        {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}};
    showMatrix(matB, "B");
    std::cout << "rank(B) = " << rank(matB) << "\n";
    std::cout << "det(B) = " << det(matB) << "\n";
    // showMatrix(inv(matB), "inv(B)"); // Error: singular matrix
    // showMatrix(pinv(matB), "pinv(B)"); // Error: not a full-rank matrix
    showMatrix(pinv2(matB), "pinv2(B)");
    showMatrix(leftDiv(matB, matb), "B \\ b"); // Warning: singular matrix
    showMatrix(matMul(pinv2(matB), matb), "pinv2(B) * b (== B \\ b)"); // Warning: singular matrix
    std::cout << "\n\n";

    i_complex_matrix matC = {
        {{1.0, 0.0}, {0.0, 1.0}},
        {{0.0, 5.0}, {1.0, 0.0}},
        {{8.0, 0.0}, {1.0, 1.0}}};
    showMatrix(matC, "C");
    std::cout << "rank(C) = " << rank(matC) << "\n";
    // std::cout << "det(C) = " << det(matC) << "\n"; // Error: not a square matrix
    // showMatrix(inv(matC), "inv(C)"); // Error: not a square matrix
    showMatrix(pinv(matC), "pinv(C)");
    showMatrix(pinv2(matC), "pinv2(C)");
    showMatrix(leftDiv(matC, matb), "C \\ b");
    showMatrix(matMul(pinv(matC), matb), "pinv(C) * b");
    showMatrix(matMul(pinv2(matC), matb), "pinv2(C) * b (== C \\ b)");
    std::cout << "\n\n";
}

int main(int argc, char **argv)
{
    realMatrixExamples();
    complexMatrixExamples();
    std::cin.get();
    return 0;
}
