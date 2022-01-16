#include "matBasic_real.hpp"
#include "matBasic_testUtil.hpp"

i_real_matrix genTestMatrixA(const std::size_t nAnt)
{
    std::size_t nrows = (nAnt - 1) * nAnt + 1;
    i_real_matrix resMat{initRealMatrix(nrows, nAnt)};
    for (std::size_t row{0}; row < nAnt; ++row)
    {
        for (std::size_t col{0}; col < nAnt; ++col)
        {
            if (col >= row)
            {
                continue;
            }
            i_float_t lambda{static_cast<i_float_t>((row + 1) * 100 + col + 1)};
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

i_real_matrix genTestMatrixA(const std::size_t nAnt, const std::size_t nEq)
{
    i_real_matrix resMat{initRealMatrix(nEq + 1, nAnt)};

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
        i_float_t lambda{static_cast<i_float_t>((row + 1) * 100 + col + 1)};
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

i_real_matrix genTestMatrixb(const std::size_t nAnt)
{
    std::size_t nrows = (nAnt - 1) * nAnt + 1;
    i_real_matrix resMat{initRealMatrix(nrows, 1)};
    resMat[nrows - 1][0] = 1.0;
    return resMat;
}

i_real_matrix genTestMatrixb(const std::size_t nAnt, const std::size_t nEq)
{
    std::size_t nrows = nEq + 1;
    i_real_matrix resMat{initRealMatrix(nrows, 1)};
    resMat[nEq][0] = 1.0;
    return resMat;
}


bool CheckIdentity(const i_real_matrix& mat)
{
    const float EPS = 1.0e-20;

    double sumErrSqr = 0;
    for (std::size_t r = 0; r < mat.size(); r++)
    {
        const i_real_vector& row = mat[r];
        for (std::size_t c = 0; c < row.size(); c++)
        {
            i_float_t diff = c == r ? row[c] - i_float_t(1.0) : row[c];
            sumErrSqr = diff * diff;
        }
    }
    
    return sumErrSqr < EPS;
}



void pinvTest(bool doLargeMatTest = true)
{
    std::cout << "\n\n******************** pinv test ********************\n\n";

    i_real_matrix matA = {
        {0.0}};
    std::cout << "rank(matA) = " << rank(matA) << "\n";
    showMatrix(matA, "matA");
    showMatrix(pinv2(matA), "pinv2(matA)");
    std::cout << "\n\n";

    i_real_matrix matB = {
        {5.0}};
    std::cout << "rank(matB) = " << rank(matB) << "\n";
    showMatrix(matB, "matB");
    showMatrix(inv(matB), "inv(matB)");
    showMatrix(pinv(matB), "pinv(matB)");
    showMatrix(pinv2(matB), "pinv2(matB)");
    std::cout << "\n\n";

    i_real_matrix matC = {
        {1.0, 2.0, 3.0, 5.0},
        {2.0, 5.0, 3.0, 6.0},
        {0.0, 4.0, 1.0, 5.0},
        {-6.0, 3.0, 9.0, 1.0}};
    std::cout << "rank(matC) = " << rank(matC) << "\n";
    showMatrix(matC, "matC");

    i_real_matrix invMatC = inv(matC);
    showMatrix(invMatC, "inv(matC)");
    bool isIdentity = CheckIdentity(matMul(matC, invMatC));
    if (!isIdentity)
        std::cout << "not identity";

    showMatrix(pinv(matC), "pinv(matC)");
    showMatrix(pinv2(matC), "pinv2(matC)");
    std::cout << "\n\n";

    i_real_matrix matD = {
        {1.0, 2.0, 3.0, 5.0},
        {2.0, 5.0, 3.0, 6.0},
        {0.0, 4.0, 1.0, 5.0},
        {-6.0, 3.0, 9.0, 1.0},
        {1.0, 1.0, 4.0, 5.0}};
    std::cout << "rank(matD) = " << rank(matD) << "\n";
    showMatrix(matD, "matD");
    showMatrix(pinv(matD), "pinv(matD)");
    showMatrix(pinv2(matD), "pinv2(matD)");
    std::cout << "\n\n";

    i_real_matrix matE = {
        {1.0, 2.0, 3.0, 5.0},
        {2.0, 5.0, 3.0, 6.0},
        {0.0, 4.0, 1.0, 5.0},
        {0.0, 4.0, 1.0, 5.0}};
    std::cout << "rank(matE) = " << rank(matE) << "\n";
    showMatrix(matE, "matE (singular)");
    showMatrix(pinv2(matE), "pinv2(matE)");
    std::cout << "\n\n";

    i_real_matrix matF = {
        {1.0, 2.0, 3.0, 0.0},
        {2.0, 5.0, 3.0, 0.0},
        {0.0, 4.0, 1.0, 0.0},
        {-6.0, 3.0, 9.0, 0.0},
        {1.0, 1.0, 4.0, 0.0}};
    std::cout << "rank(matF) = " << rank(matF) << "\n";
    showMatrix(matF, "matF (singular)");
    showMatrix(pinv2(matF), "pinv2(matF)");
    std::cout << "\n\n";

    i_real_matrix matG = {
        {0.0, 2.0, 3.0, 5.0},
        {2.0, 0.0, 3.0, 6.0},
        {0.0, 4.0, 0.0, 5.0},
        {-6.0, 3.0, 9.0, 0.0}};
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
    // i_real_matrix largeA = genTestMatrixA(nAnt, nEq);
    // i_real_matrix largeb = genTestMatrixb(nAnt, nEq);
    i_real_matrix largeA = genTestMatrixA(nAnt);
    i_real_matrix largeb = genTestMatrixb(nAnt);
    i_real_matrix largeA_pinv;
    i_real_matrix largex;

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
    i_real_matrix matA = {
        {5.0}};
    showMatrix(matA, "matA", true);
    std::cout << "det(matA) = " << det(matA) << "\n\n\n";

    i_real_matrix matB = {
        {0.0, 1.0},
        {4.0, 0.0}};
    showMatrix(matB, "matB", true);
    std::cout << "det(matB) = " << det(matB) << "\n\n\n";

    i_real_matrix matC = {
        {0.0, 1.0, 1.0},
        {4.0, 0.0, 1.0},
        {8.0, 1.0, 0.0}};
    showMatrix(matC, "matC", true);
    std::cout << "det(matC) = " << det(matC) << "\n\n\n";

    i_real_matrix matD = {
        {0.0, 1.0, 1.0, 4.0},
        {4.0, 0.0, 1.0, 9.0},
        {8.0, 1.0, 0.0, 0.0},
        {1.0, 1.0, 1.0, 0.0}};
    showMatrix(matD, "matD", true);
    std::cout << "det(matD) = " << det(matD) << "\n\n\n";

    i_real_matrix matE = {
        {0.0, 1.0, 1.0, 4.0, 5.0},
        {4.0, 0.0, 1.0, 9.0, 1.0},
        {8.0, 1.0, 0.0, 0.0, 9.0},
        {1.0, 1.0, 1.0, 0.0, 4.0},
        {1.0, 4.0, 8.0, 1.0, 0.0}};
    showMatrix(matE, "matE", true);
    std::cout << "det(matE) = " << det(matE) << "\n\n\n";

    i_real_matrix matF = {
        {0.0, 1.0, 1.0, 4.0, 5.0, 1.0},
        {4.0, 0.0, 1.0, 9.0, 1.0, 9.0},
        {8.0, 1.0, 0.0, 0.0, 9.0, 3.0},
        {1.0, 1.0, 1.0, 0.0, 4.0, 5.0},
        {1.0, 4.0, 8.0, 1.0, 0.0, 0.0},
        {9.0, 3.0, 1.0, 1.0, 4.0, 0.0}};
    showMatrix(matF, "matF", true);
    std::cout << "det(matF) = " << det(matF) << "\n\n\n";

    // std::size_t nAnt = 64;
    // std::size_t nEq = 960;
    // i_real_matrix largeA = genTestMatrixA(nAnt, nEq);
    // i_real_matrix largeG = matMul(transpose(largeA), largeA);
    // TestTimer timer;
    // timer.tic();
    // std::cout << "det(largeG) = " << det(largeG) << "\n";
    // timer.toc("large matrix determinant");
}

int main(int argc, char **argv)
{
    pinvTest(true);
    determinantTest();
    std::cin.get();
    return 0;
}
