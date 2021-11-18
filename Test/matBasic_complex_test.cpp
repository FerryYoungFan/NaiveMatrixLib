#include "matBasic_complex.hpp"
#include "matBasic_testUtil.hpp"

i_complex_matrix genTestMatrixA(const std::size_t nAnt)
{
    std::size_t nrows = (nAnt - 1) * nAnt + 1;
    i_complex_matrix resMat{initComplexMatrix(nrows, nAnt)};
    for (std::size_t row{0}; row < nAnt; ++row)
    {
        for (std::size_t col{0}; col < nAnt; ++col)
        {
            if (col >= row)
            {
                continue;
            }
            // i_complex_t lambda{static_cast<i_float_t>((row + 1) * 10 + col + 1)};
            i_complex_t lambda{static_cast<i_float_t>(row + 1), static_cast<i_float_t>(col + 1)};
            std::size_t row_this = col * (nAnt - 1) + row - (row > col);
            std::size_t row_that = row * (nAnt - 1) + col;
            resMat[row_this][row] = lambda;
            resMat[row_this][col] = -1.0;
            resMat[row_that][col] = lambda;
            resMat[row_that][row] = -1.0;
        }
    }
    resMat[nrows - 1][0] = 1.0;
    return resMat;
}

i_complex_matrix genTestMatrixA(const std::size_t nAnt, const std::size_t nEq)
{
    i_complex_matrix resMat{initComplexMatrix(nEq + 1, nAnt)};

    if (nAnt > nEq)
    {
        std::cout << "genTestMatrixA Error: nEq should >= nAnt\n";
        return resMat;
    }

    std::size_t maxRow = (nAnt - 1) * nAnt;
    if (nEq > maxRow)
    {
        std::cout << "genTestMatrixA Error: nEq should <= " << maxRow << "\n";
        return resMat;
    }

    std::size_t antCount{0};
    std::size_t col{0};
    std::size_t roundCount{0};
    for (std::size_t row{0}; row < nEq; ++row)
    {
        resMat[row][antCount] = -1.0;
        col = (antCount + roundCount + (antCount >= roundCount)) % nAnt;
        //  lambda{static_cast<i_float_t>((col + 1)*10 + antCount + 1), 0.0};
        i_complex_t lambda{static_cast<i_float_t>(col + 1), static_cast<i_float_t>(antCount + 1)};
        resMat[row][col] = lambda;
        ++antCount;
        if (antCount == nAnt)
        {
            antCount = 0;
            ++roundCount;
        }
    }
    resMat[nEq][0] = 1.0;
    return resMat;
}

i_complex_matrix genTestMatrixb(const std::size_t nAnt)
{
    std::size_t nrows = (nAnt - 1) * nAnt + 1;
    i_complex_matrix resMat{initComplexMatrix(nrows, 1)};
    resMat[nrows - 1][0] = 1.0;
    return resMat;
}

i_complex_matrix genTestMatrixb(const std::size_t nAnt, const std::size_t nEq)
{
    std::size_t nrows = nEq + 1;
    i_complex_matrix resMat{initComplexMatrix(nrows, 1)};
    resMat[nEq][0] = 1.0;
    return resMat;
}

void pinvTest(bool doLargeMatTest = true)
{
    std::cout << "\n\n******************** pinv test ********************\n\n";

    i_complex_matrix matA = {
        {{0.0, 0.0}}};
    std::cout << "rank(matA) = " << rank(matA) << "\n";
    showMatrix(matA, "matA");
    showMatrix(pinv2(matA), "pinv2(matA)");
    std::cout << "\n\n";

    i_complex_matrix matB = {
        {{5.0, 3.0}}};
    std::cout << "rank(matB) = " << rank(matB) << "\n";
    showMatrix(matB, "matB");
    showMatrix(inv(matB), "inv(matB)");
    showMatrix(pinv(matB), "pinv(matB)");
    showMatrix(pinv2(matB), "pinv2(matB)");
    std::cout << "\n\n";

    i_complex_matrix matC = {
        {{1.0, 0.0}, {2.0, 0.0}, {3.0, 2.0}, {5.0, 0.0}},
        {{2.0, 0.0}, {5.0, 0.0}, {3.0, 0.0}, {6.0, 0.0}},
        {{0.0, 7.0}, {4.0, 0.0}, {3.0, 0.0}, {1.0, 0.0}},
        {{-6.0, 0.0}, {3.0, 0.0}, {1.0, 0.0}, {9.0, 1.0}}};
    std::cout << "rank(matC) = " << rank(matC) << "\n";
    showMatrix(matC, "matC");
    showMatrix(inv(matC), "inv(matC)");
    showMatrix(pinv(matC), "pinv(matC)");
    showMatrix(pinv2(matC), "pinv2(matC)");
    std::cout << "\n\n";

    i_complex_matrix matD = {
        {{1.0, 0.0}, {2.0, 0.0}, {3.0, 2.0}, {5.0, 0.0}},
        {{2.0, 0.0}, {5.0, 0.0}, {3.0, 0.0}, {6.0, 0.0}},
        {{0.0, 7.0}, {4.0, 0.0}, {3.0, 0.0}, {1.0, 0.0}},
        {{-6.0, 0.0}, {3.0, 0.0}, {1.0, 0.0}, {9.0, 1.0}},
        {{3.0, 0.0}, {2.0, 0.0}, {1.0, 5.0}, {0.0, 1.0}}};
    std::cout << "rank(matD) = " << rank(matD) << "\n";
    showMatrix(matD, "matD");
    showMatrix(pinv(matD), "pinv(matD)");
    showMatrix(pinv2(matD), "pinv2(matD)");
    std::cout << "\n\n";

    i_complex_matrix matE = {
        {{1.0, 0.0}, {2.0, 0.0}, {3.0, 2.0}, {5.0, 0.0}},
        {{2.0, 0.0}, {5.0, 0.0}, {3.0, 0.0}, {6.0, 0.0}},
        {{0.0, 7.0}, {4.0, 0.0}, {3.0, 0.0}, {1.0, 0.0}},
        {{0.0, 7.0}, {4.0, 0.0}, {3.0, 0.0}, {1.0, 0.0}}};
    std::cout << "rank(matE) = " << rank(matE) << "\n";
    showMatrix(matE, "matE (singular)");
    showMatrix(pinv2(matE), "pinv2(matE)");
    std::cout << "\n\n";

    i_complex_matrix matF = {
        {{1.0, 0.0}, {2.0, 0.0}, {3.0, 2.0}, {0.0, 0.0}},
        {{2.0, 0.0}, {5.0, 0.0}, {3.0, 0.0}, {0.0, 0.0}},
        {{0.0, 7.0}, {4.0, 0.0}, {3.0, 0.0}, {0.0, 0.0}},
        {{-6.0, 0.0}, {3.0, 0.0}, {1.0, 0.0}, {0.0, 0.0}},
        {{3.0, 0.0}, {2.0, 0.0}, {1.0, 5.0}, {0.0, 0.0}}};
    std::cout << "rank(matF) = " << rank(matF) << "\n";
    showMatrix(matF, "matF (singular)");
    showMatrix(pinv2(matF), "pinv2(matF)");
    std::cout << "\n\n";

    i_complex_matrix matG = {
        {{0.0, 0.0}, {2.0, 0.0}, {3.0, 2.0}, {5.0, 0.0}},
        {{2.0, 0.0}, {0.0, 0.0}, {3.0, 0.0}, {6.0, 0.0}},
        {{0.0, 7.0}, {4.0, 0.0}, {0.0, 0.0}, {1.0, 0.0}},
        {{-6.0, 0.0}, {3.0, 0.0}, {1.0, 0.0}, {0.0, 0.0}}};
    std::cout << "rank(matG) = " << rank(matG) << "\n";
    showMatrix(matG, "matG (diag zero)");
    showMatrix(inv(matG), "inv(matG)");
    showMatrix(pinv(matG), "pinv(matG)");
    showMatrix(pinv2(matG), "pinv2(matG)");
    std::cout << "\n\n";

    if (!doLargeMatTest)
    {
        return;
    }

    std::size_t nAnt = 64;
    std::size_t nEq = 960;
    i_complex_matrix largeA = genTestMatrixA(nAnt, nEq);
    i_complex_matrix largeb = genTestMatrixb(nAnt, nEq);
    // i_complex_matrix largeA = genTestMatrixA(nAnt);
    // i_complex_matrix largeb = genTestMatrixb(nAnt);
    i_complex_matrix largeA_pinv;
    i_complex_matrix largex;

    TestTimer timer;
    timer.tic();
    largeA_pinv = pinv(largeA);
    largex = matMul(largeA_pinv, largeb);
    showMatrix(largex, "mat x");
    timer.toc("pinv1 method");

    timer.tic();
    largex = leftDiv(largeA, largeb);
    showMatrix(largex, "mat x");
    timer.toc("pinv2 method");
}

void determinantTest()
{
    std::cout << "\n\n******************** determinant test ********************\n\n";
    i_complex_matrix matA = {
        {{5.0, 3.0}}};
    showMatrix(matA, "matA", true);
    std::cout << "det(matA) = " << det(matA) << "\n\n\n";

    i_complex_matrix matB = {
        {{1.0, 0.0}, {2.0, 0.0}},
        {{2.0, 0.0}, {0.0, 5.0}}};
    showMatrix(matB, "matB", true);
    std::cout << "det(matB) = " << det(matB) << "\n\n\n";

    i_complex_matrix matC = {
        {{1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}},
        {{2.0, 0.0}, {5.0, 0.0}, {3.0, 0.0}},
        {{0.0, 7.0}, {4.0, 0.0}, {3.0, 0.0}}};
    showMatrix(matC, "matC", true);
    std::cout << "det(matC) = " << det(matC) << "\n\n\n";

    i_complex_matrix matD = {
        {{1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {5.0, 0.0}},
        {{2.0, 0.0}, {5.0, 0.0}, {3.0, 0.0}, {6.0, 0.0}},
        {{0.0, 7.0}, {4.0, 0.0}, {3.0, 0.0}, {1.0, 0.0}},
        {{-6.0, 0.0}, {3.0, 0.0}, {1.0, 0.0}, {9.0, 0.0}}};
    showMatrix(matD, "matD", true);
    std::cout << "det(matD) = " << det(matD) << "\n\n\n";

    i_complex_matrix matE = {
        {{1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {5.0, 0.0}, {1.0, 0.0}},
        {{2.0, 0.0}, {5.0, 0.0}, {3.0, 0.0}, {6.0, 0.0}, {3.0, 0.0}},
        {{0.0, 7.0}, {4.0, 0.0}, {3.0, 0.0}, {1.0, 0.0}, {1.0, 0.0}},
        {{-6.0, 0.0}, {3.0, 0.0}, {1.0, 0.0}, {9.0, 0.0}, {4.0, 0.0}},
        {{-3.0, 0.0}, {3.0, 5.0}, {7.0, 0.0}, {0.0, 0.0}, {1.0, 0.0}}};
    showMatrix(matE, "matE", true);
    std::cout << "det(matE) = " << det(matE) << "\n\n\n";

    i_complex_matrix matF = {
        {{9.0, 12.0}, {2.0, 0.0}, {3.0, 0.0}, {5.0, 0.0}, {1.0, 0.0}, {1.0, 7.0}},
        {{2.0, 0.0}, {5.0, 0.0}, {3.0, 0.0}, {6.0, 0.0}, {3.0, 0.0}, {0.0, 0.0}},
        {{0.0, 7.0}, {4.0, 0.0}, {3.0, 0.0}, {1.0, 0.0}, {1.0, 0.0}, {0.0, 0.0}},
        {{-6.0, 0.0}, {3.0, 0.0}, {1.0, 0.0}, {9.0, 0.0}, {4.0, 0.0}, {5.0, 0.0}},
        {{-3.0, 0.0}, {3.0, 5.0}, {7.0, 0.0}, {0.0, 0.0}, {1.0, 0.0}, {1.0, 0.0}},
        {{1.0, 0.0}, {5.0, 5.0}, {1.0, 0.0}, {0.0, 1.0}, {0.0, 0.0}, {9.0, 0.0}}};
    showMatrix(matF, "matF", true);
    std::cout << "det(matF) = " << det(matF) << "\n\n\n";

    std::size_t nAnt = 64;
    std::size_t nEq = 960;
    i_complex_matrix largeA = genTestMatrixA(nAnt, nEq);
    i_complex_matrix largeG = matMul(transpose(largeA), largeA);
    TestTimer timer;
    timer.tic();
    std::cout << "det(largeG) = " << det(largeG) << "\n";
    timer.toc("large matrix determinant");
}

int main(int argc, char **argv)
{
    pinvTest(true);
    determinantTest();
    std::cin.get();
    return 0;
}
