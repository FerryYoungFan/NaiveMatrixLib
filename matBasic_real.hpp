#pragma once
/*
Mini Matrix Tools for Real Matrix
Vesion: 0.9 - 2021.11.18 by Fanseline
https://ferryyoungfan.github.io

Main Function List:
[1] rank:    Matrix rank (Cholesky decomposition)
[2] det:     Matrix determinant calculation
[3] inv:     LU decomposition-based matrix inversion
[4] pinv:    pinv(G) = inv(G' * G) * G' (WARNING: full-rank matrix only!)
[5] pinv2:   Moore-Penrose pseudoinversion (same as pinv(G) in MATLAB)
[6] leftDiv: x = A \ b, using Moore-Penrose pinv, NOT same as MATLAB for a singular matrix

Reference:
[*1] Pierre Courrieu, Fast Computation of Moore-Penrose Inverse Matrices, https://arxiv.org/abs/0804.4809
[*2] Permute Sign Calculation, page5 https://www.math.rutgers.edu/docman-lister/math-main/academics/course-materials/250/assignments/1493-250c-lab3-sakai-pdf/file
[*3] LU Decomposition C++ Implementation, https://blog.csdn.net/xx_123_1_rj/article/details/39553809
[*4] LU Decomposition, https://www.math.ucdavis.edu/~linear/old/notes11.pdf
*/
#include <iostream>
#include <cmath>
#include <vector>

using i_float_t = double; // using i_float_t = float; // Notice: Do NOT use int type!
using i_real_vector = std::vector<i_float_t>;
using i_real_matrix = std::vector<i_real_vector>;

// Simply print real matrix with description, can be either block or MATLAB format.
void showMatrix(const i_real_matrix &matG, const char *describe = nullptr, bool matlabFormat = false)
{
    const std::size_t nrows{matG.size()}, ncols{matG[0].size()};
    matlabFormat = describe && matlabFormat;
    if (describe)
    {
        matlabFormat ? std::cout << describe << " = [" : std::cout << describe << " : " << nrows << " x " << ncols << " Real Matrix:\n";
    }
    for (std::size_t row{0}; row < nrows; ++row)
    {
        if (!matlabFormat)
        {
            std::cout << "    row[" << row + 1 << "]: ";
        }
        for (std::size_t col{0}; col < ncols; ++col)
        {
            std::cout << matG[row][col];
            if (col + 1 < ncols)
            {
                std::cout << ",  ";
            }
            else
            {
                matlabFormat ? std::cout << "; " : std::cout << ";\n";
            }
        }
    }
    matlabFormat ? std::cout << "];\n" : std::cout << "\n";
}

// Generate and fill an nrows x ncols matrix, fill with given value (zero by default)
i_real_matrix initRealMatrix(const std::size_t nrows, const std::size_t ncols, const i_float_t initValue = 0.0)
{
    return i_real_matrix(nrows, i_real_vector(ncols, initValue));
}

// Matrix transpose
i_real_matrix transpose(const i_real_matrix &matG)
{
    const std::size_t nrows{matG.size()}, ncols{matG[0].size()};
    i_real_matrix matGt = initRealMatrix(ncols, nrows);
    std::size_t i{0}, j{0};
    for (i = 0; i < nrows; ++i)
    {
        for (j = 0; j < ncols; ++j)
        {
            matGt[j][i] = matG[i][j];
        }
    }
    return matGt;
}

// Matrix multiplication O(n^3) naive implementation
i_real_matrix matMul(const i_real_matrix &matA, const i_real_matrix &matB)
{
    const std::size_t nrowsA{matA.size()}, ncolsA{matA[0].size()}, nrowsB{matB.size()}, ncolsB{matB[0].size()};
    i_real_matrix resMat;
    if (ncolsA != nrowsB)
    {
        std::cout << "Error when using matMul: dimension not match.\n";
        return resMat;
    }
    resMat = initRealMatrix(nrowsA, ncolsB);
    std::size_t i{0}, j{0}, k{0};
    for (i = 0; i < nrowsA; ++i)
    {
        for (j = 0; j < ncolsB; ++j)
        {
            for (k = 0; k < ncolsA; ++k)
            {
                resMat[i][j] += matA[i][k] * matB[k][j];
            }
        }
    }
    return resMat;
}

// Calculate matrix rank (Cholesky decomposition) [*1]
std::size_t rank(const i_real_matrix &matG, const i_float_t tolerance = 1.0e-9)
{
    const std::size_t nrows{matG.size()}, ncols{matG[0].size()};
    std::size_t nSize{ncols};
    std::size_t i{0}, j{0}, k{0};

    i_real_matrix matA;
    if (nrows < nSize)
    {
        // A = G * G'
        nSize = nrows;
        matA = initRealMatrix(nSize, nSize);
        for (i = 0; i < nSize; ++i)
        {
            for (j = 0; j < nSize; ++j)
            {
                for (k = 0; k < ncols; ++k)
                {
                    matA[i][j] += matG[i][k] * matG[j][k];
                }
            }
        }
    }
    else
    {
        // A = G' * G
        matA = initRealMatrix(nSize, nSize);
        for (i = 0; i < nSize; ++i)
        {
            for (j = 0; j < nSize; ++j)
            {
                for (k = 0; k < ncols; ++k)
                {
                    matA[i][j] += matG[k][i] * matG[k][j];
                }
            }
        }
    }

    // Full rank Cholesky decomposition of A
    i_float_t tol{std::abs(matA[0][0])};
    for (i = 0; i < nSize; ++i)
    {
        if (matA[i][i] > 0)
        {
            const i_float_t temp{std::abs(matA[i][i])};
            if (temp < tol)
            {
                tol = temp;
            }
        }
    }
    tol *= tolerance;

    i_real_matrix matL = initRealMatrix(nSize, nSize);
    std::size_t rankA{0};
    for (k = 0; k < nSize; ++k)
    {
        for (i = k; i < nSize; ++i)
        {
            matL[i][rankA] = matA[i][k];
            for (j = 0; j < rankA; ++j)
            {
                matL[i][rankA] -= matL[i][j] * matL[k][j];
            }
        }
        if (matL[k][rankA] > tol)
        {
            matL[k][rankA] = std::sqrt(matL[k][rankA]);
            if (k < nSize)
            {
                for (j = k + 1; j < nSize; ++j)
                {
                    matL[j][rankA] /= matL[k][rankA];
                }
            }
            ++rankA;
        }
    }
    return rankA; // rank(G) = rank(A)
}

// LU decomposition-based matrix determinant calculation [*2][*3][*4]
i_float_t det(const i_real_matrix &matG)
{
    const std::size_t nrows{matG.size()}, ncols{matG[0].size()};
    i_float_t detG = 0.0;
    if (nrows != ncols)
    {
        std::cout << "Error when using det: matrix is not square.\n";
        return detG;
    }

    const std::size_t nSize{nrows};
    std::size_t i{0}, j{0}, k{0};

    // ******************** Step 1: row permutation (swap diagonal zeros) ********************
    i_real_matrix matLU;
    std::vector<std::size_t> permuteLU; // Permute vector
    bool changeSign{false};

    for (i = 0; i < nSize; ++i)
    {
        permuteLU.push_back(i);
    }

    for (j = 0; j < nSize; ++j)
    {
        i_float_t maxv{0.0};
        for (i = j; i < nSize; ++i)
        {
            const i_float_t currentv{std::abs(matG[permuteLU[i]][j])};
            if (currentv > maxv)
            {
                maxv = currentv;
                if (permuteLU[i] != permuteLU[j]) // swap rows
                {
                    changeSign = !changeSign;
                    const std::size_t tmp{permuteLU[j]};
                    permuteLU[j] = permuteLU[i];
                    permuteLU[i] = tmp;
                }
            }
        }
    }

    for (i = 0; i < nSize; ++i)
    {
        matLU.push_back(matG[permuteLU[i]]);
    }

    // ******************** Step 2: LU decomposition (save both L & U in matLU) ********************
    if (matLU[0][0] == 0.0)
    {
        return detG; // Singular matrix, det(G) = 0
    }

    for (i = 1; i < nSize; ++i)
    {
        matLU[i][0] /= matLU[0][0];
    }

    for (i = 1; i < nSize; ++i)
    {
        for (j = i; j < nSize; ++j)
        {
            for (k = 0; k < i; ++k)
            {
                matLU[i][j] -= matLU[i][k] * matLU[k][j]; // Calculate U matrix
            }
        }
        if (matLU[i][i] == 0.0)
        {
            return detG; // Singular matrix, det(G) = 0
        }
        for (k = i + 1; k < nSize; ++k)
        {
            for (j = 0; j < i; ++j)
            {
                matLU[k][i] -= matLU[k][j] * matLU[j][i]; // Calculate L matrix
            }
            matLU[k][i] /= matLU[i][i];
        }
    }

    detG = 1.0;
    if (changeSign)
    {
        detG = -1.0; // Change the sign of the determinant
    }
    for (i = 0; i < nSize; ++i)
    {
        detG *= matLU[i][i]; // det(G) = det(L) * det(U). For triangular matrices, det(L) = prod(diag(L)) = 1, det(U) = prod(diag(U)), so det(G) = prod(diag(U))
    }

    return detG;
}

// LU decomposition-based matrix inversion [*3][*4]
i_real_matrix inv(const i_real_matrix &matG, const bool usePermute = true)
{
    const std::size_t nrows{matG.size()}, ncols{matG[0].size()};
    i_real_matrix matLU;
    if (nrows != ncols)
    {
        std::cout << "Error when using inv: matrix is not square.\n";
        return matLU;
    }

    const std::size_t nSize{nrows};
    std::size_t i{0}, j{0}, k{0};

    // ******************** Step 1: row permutation (swap diagonal zeros) ********************
    std::vector<std::size_t> permuteLU; // Permute vector
    for (i = 0; i < nSize; ++i)
    {
        permuteLU.push_back(i); // Push back row index
    }

    if (usePermute) // Sort rows by pivot element
    {
        for (j = 0; j < nSize; ++j)
        {
            i_float_t maxv{0.0};
            for (i = j; i < nSize; ++i)
            {
                const i_float_t currentv{std::abs(matG[permuteLU[i]][j])};
                if (currentv > maxv) // Swap rows
                {
                    maxv = currentv;
                    const std::size_t tmp{permuteLU[j]};
                    permuteLU[j] = permuteLU[i];
                    permuteLU[i] = tmp;
                }
            }
        }
        for (i = 0; i < nSize; ++i)
        {
            matLU.push_back(matG[permuteLU[i]]); // Make a permuted matrix with new row order
        }
    }
    else
    {
        matLU = i_real_matrix(matG); // Simply duplicate matrix
    }

    // ******************** Step 2: LU decomposition (save both L & U in matLU) ********************
    if (matLU[0][0] == 0.0)
    {
        std::cout << "Warning when using inv: matrix is singular.\n";
        matLU.clear();
        return matLU;
    }
    for (i = 1; i < nSize; ++i)
    {
        matLU[i][0] /= matLU[0][0]; // Initialize first column of L matrix
    }
    for (i = 1; i < nSize; ++i)
    {
        for (j = i; j < nSize; ++j)
        {
            for (k = 0; k < i; ++k)
            {
                matLU[i][j] -= matLU[i][k] * matLU[k][j]; // Calculate U matrix
            }
        }
        if (matLU[i][i] == 0.0)
        {
            std::cout << "Warning when using inv: matrix is singular.\n";
            matLU.clear();
            return matLU;
        }
        for (k = i + 1; k < nSize; ++k)
        {
            for (j = 0; j < i; ++j)
            {
                matLU[k][i] -= matLU[k][j] * matLU[j][i]; // Calculate L matrix
            }
            matLU[k][i] /= matLU[i][i];
        }
    }

    // ******************** Step 3: L & U inversion (save both L^-1 & U^-1 in matLU_inv) ********************
    i_real_matrix matLU_inv = initRealMatrix(nSize, nSize);

    // matL inverse & matU inverse
    for (i = 0; i < nSize; ++i)
    {
        // L matrix inverse, omit diagonal ones
        matLU_inv[i][i] = 1.0;
        for (k = i + 1; k < nSize; ++k)
        {
            for (j = i; j <= k - 1; ++j)
            {
                matLU_inv[k][i] -= matLU[k][j] * matLU_inv[j][i];
            }
        }
        // U matrix inverse
        matLU_inv[i][i] = 1.0 / matLU[i][i];
        for (k = i; k > 0; --k)
        {
            for (j = k; j <= i; ++j)
            {
                matLU_inv[k - 1][i] -= matLU[k - 1][j] * matLU_inv[j][i];
            }
            matLU_inv[k - 1][i] /= matLU[k - 1][k - 1];
        }
    }

    // ******************** Step 4: Calculate G^-1 = U^-1 * L^-1 ********************
    // Lower part product
    for (i = 1; i < nSize; ++i)
    {
        for (j = 0; j < i; ++j)
        {
            const std::size_t jp{permuteLU[j]}; // Permute column back
            matLU[i][jp] = 0.0;
            for (k = i; k < nSize; ++k)
            {
                matLU[i][jp] += matLU_inv[i][k] * matLU_inv[k][j];
            }
        }
    }
    // Upper part product
    for (i = 0; i < nSize; ++i)
    {
        for (j = i; j < nSize; ++j)
        {
            const std::size_t jp{permuteLU[j]}; // Permute column back
            matLU[i][jp] = matLU_inv[i][j];
            for (k = j + 1; k < nSize; ++k)
            {
                matLU[i][jp] += matLU_inv[i][k] * matLU_inv[k][j];
            }
        }
    }
    return matLU; // Reused matLU as a result container
}

// Classic pseudoinversion pinv(G) = inv(G' * G) * G' (WARNING: full-rank matrix only!)
i_real_matrix pinv(const i_real_matrix &matG)
{
    i_real_matrix matGt = transpose(matG);
    i_real_matrix matGtG_inv = inv(matMul(matGt, matG));
    return matMul(matGtG_inv, matGt);
}

// Moore-Penrose pseudoinversion (same as pinv(G) in MATLAB) [*1]
i_real_matrix pinv2(const i_real_matrix &matG, const i_float_t tolerance = 1.0e-9)
{
    bool useTranspose{false};
    const std::size_t nrows{matG.size()}, ncols{matG[0].size()};
    std::size_t nSize{ncols};

    i_real_matrix matA, matGt;
    matGt = transpose(matG);
    if (nrows < nSize)
    {
        useTranspose = true;
        nSize = nrows;
        matA = matMul(matG, matGt); // A = G * G'
    }
    else
    {
        matA = matMul(matGt, matG); // A = G' * G
    }

    // Full rank Cholesky decomposition of A
    std::size_t i{0}, j{0}, k{0};

    i_float_t tol{std::abs(matA[0][0])};
    for (i = 0; i < nSize; ++i)
    {
        if (matA[i][i] > 0)
        {
            const i_float_t temp{matA[i][i]};
            if (temp < tol)
            {
                tol = temp;
            }
        }
    }
    tol *= tolerance;

    i_real_matrix matL = initRealMatrix(nSize, nSize);
    std::size_t rankA{0};
    for (k = 0; k < nSize; ++k)
    {
        for (i = k; i < nSize; ++i)
        {
            matL[i][rankA] = matA[i][k];
            for (j = 0; j < rankA; ++j)
            {
                matL[i][rankA] -= matL[i][j] * matL[k][j];
            }
        }
        if (matL[k][rankA] > tol)
        {
            matL[k][rankA] = std::sqrt(matL[k][rankA]);
            if (k < nSize)
            {
                for (j = k + 1; j < nSize; ++j)
                {
                    matL[j][rankA] /= matL[k][rankA];
                }
            }
            ++rankA;
        }
    }

    if (rankA == 0)
    {
        return matGt; // All-zero matrix's transpose
    }

    // Slice L = L(:, 0:r);
    for (i = 0; i < nSize; ++i)
    {
        for (k = 0; k < nSize - rankA; ++k)
        {
            matL[i].pop_back();
        }
    }

    // Generalized inverse
    i_real_matrix matLt = transpose(matL);
    i_real_matrix matM = inv(matMul(matLt, matL), false);   // M = inv(L' * L)
    matA = matMul(matMul(matMul(matL, matM), matM), matLt); // A = L * M * M * L'

    if (useTranspose)
    {
        return matMul(matGt, matA); // pinv(G) = G' * (L * M * M * L')
    }
    return matMul(matA, matGt); // pinv(G) = (L * M * M * L') * G'
}

// Calculate left division x = A \ b, using Moore-Penrose pinv, NOT same as MATLAB for a singular matrix
i_real_matrix leftDiv(const i_real_matrix &matA, const i_real_matrix &matb)
{
    i_real_matrix matx;
    if (matA.size() != matb.size())
    {
        std::cout << "Error when using leftDiv: row size not match.\n";
        return matx;
    }
    matx = matMul(pinv2(matA), matb); // x = A \ b = pinv(A) * b
    return matx;
}
