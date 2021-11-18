#pragma once
/*
Mini Matrix Tools for Complex Matrix
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
#include <complex>
#include <vector>

using i_float_t = double; // using i_float_t = float; // Notice: Do NOT use int type!
using i_complex_t = std::complex<i_float_t>;
using i_complex_vector = std::vector<i_complex_t>;
using i_complex_matrix = std::vector<i_complex_vector>;

// Simply print complex matrix with description, can be either block or MATLAB format.
void showMatrix(const i_complex_matrix &matG, const char *describe = nullptr, bool matlabFormat = false)
{
    const std::size_t nrows{matG.size()}, ncols{matG[0].size()};
    matlabFormat = describe && matlabFormat;
    if (describe)
    {
        matlabFormat ? std::cout << describe << " = [" : std::cout << describe << " : " << nrows << " x " << ncols << " Complex Matrix:\n";
    }
    for (std::size_t row{0}; row < nrows; ++row)
    {
        if (!matlabFormat)
        {
            std::cout << "    row[" << row + 1 << "]: ";
        }
        for (std::size_t col{0}; col < ncols; ++col)
        {
            std::cout << matG[row][col].real();
            const i_float_t imag{matG[row][col].imag()};
            if (imag != 0)
            {
                std::cout << (imag > 0 ? "+" : "-") << std::abs(imag) << "i";
            }
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
i_complex_matrix initComplexMatrix(const std::size_t nrows, const std::size_t ncols, const i_complex_t initValue = i_complex_t{0.0, 0.0})
{
    return i_complex_matrix(nrows, i_complex_vector(ncols, initValue));
}

// Conjugate transpose (a.k.a. Hermitian transpose, G' = G^H = conj(G^T))
i_complex_matrix transpose(const i_complex_matrix &matG)
{
    const std::size_t nrows{matG.size()}, ncols{matG[0].size()};
    i_complex_matrix matGt = initComplexMatrix(ncols, nrows);
    std::size_t i{0}, j{0};
    for (i = 0; i < nrows; ++i)
    {
        for (j = 0; j < ncols; ++j)
        {
            matGt[j][i] = std::conj(matG[i][j]); // Hermitian transpose
        }
    }
    return matGt;
}

// Matrix multiplication O(n^3) naive implementation
i_complex_matrix matMul(const i_complex_matrix &matA, const i_complex_matrix &matB)
{
    const std::size_t nrowsA{matA.size()}, ncolsA{matA[0].size()}, nrowsB{matB.size()}, ncolsB{matB[0].size()};
    i_complex_matrix resMat;
    if (ncolsA != nrowsB)
    {
        std::cout << "Error when using matMul: dimension not match.\n";
        return resMat;
    }
    resMat = initComplexMatrix(nrowsA, ncolsB);
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
std::size_t rank(const i_complex_matrix &matG, const i_float_t tolerance = 1.0e-9)
{
    const std::size_t nrows{matG.size()}, ncols{matG[0].size()};
    std::size_t nSize{ncols};
    std::size_t i{0}, j{0}, k{0};

    i_complex_matrix matA;
    if (nrows < nSize)
    {
        // A = G * G'
        nSize = nrows;
        matA = initComplexMatrix(nSize, nSize);
        for (i = 0; i < nSize; ++i)
        {
            for (j = 0; j < nSize; ++j)
            {
                for (k = 0; k < ncols; ++k)
                {
                    matA[i][j] += matG[i][k] * std::conj(matG[j][k]);
                }
            }
        }
    }
    else
    {
        // A = G' * G
        matA = initComplexMatrix(nSize, nSize);
        for (i = 0; i < nSize; ++i)
        {
            for (j = 0; j < nSize; ++j)
            {
                for (k = 0; k < ncols; ++k)
                {
                    matA[i][j] += std::conj(matG[k][i]) * matG[k][j];
                }
            }
        }
    }

    // Full rank Cholesky decomposition of A
    i_float_t tol{std::abs(matA[0][0])};
    for (i = 0; i < nSize; ++i)
    {
        if (matA[i][i].real() > 0)
        {
            const i_float_t temp{std::abs(matA[i][i])};
            if (temp < tol)
            {
                tol = temp;
            }
        }
    }
    tol *= tolerance;

    i_complex_matrix matL = initComplexMatrix(nSize, nSize);
    std::size_t rankA{0};
    for (k = 0; k < nSize; ++k)
    {
        for (i = k; i < nSize; ++i)
        {
            matL[i][rankA] = matA[i][k];
            for (j = 0; j < rankA; ++j)
            {
                matL[i][rankA] -= matL[i][j] * std::conj(matL[k][j]);
            }
        }
        if (matL[k][rankA].real() > tol)
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
i_complex_t det(const i_complex_matrix &matG)
{
    const std::size_t nrows{matG.size()}, ncols{matG[0].size()};
    i_complex_t detG = 0.0;
    if (nrows != ncols)
    {
        std::cout << "Error when using det: matrix is not square.\n";
        return detG;
    }

    const std::size_t nSize{nrows};
    std::size_t i{0}, j{0}, k{0};

    // ******************** Step 1: row permutation (swap diagonal zeros) ********************
    i_complex_matrix matLU;
    std::vector<std::size_t> permuteLU; // Permute vector
    bool changeSign{false};

    for (i = 0; i < nSize; ++i)
    {
        permuteLU.push_back(i);
    }

    for (j = 0; j < nSize; ++j)
    {
        i_float_t maxv{0.0};
        i_float_t currentv{0.0};
        for (i = j; i < nSize; ++i)
        {
            if (matG[permuteLU[i]][j].real() != 0)
            {
                currentv = std::abs(matG[permuteLU[i]][j].real());
            }
            else
            {
                currentv = std::abs(matG[permuteLU[i]][j].imag());
            }
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
i_complex_matrix inv(const i_complex_matrix &matG, const bool usePermute = true)
{
    const std::size_t nrows{matG.size()}, ncols{matG[0].size()};
    i_complex_matrix matLU;
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
            i_float_t currentv{0.0};
            for (i = j; i < nSize; ++i)
            {
                if (matG[permuteLU[i]][j].real() != 0)
                {
                    currentv = std::abs(matG[permuteLU[i]][j].real());
                }
                else
                {
                    currentv = std::abs(matG[permuteLU[i]][j].imag());
                }
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
        matLU = i_complex_matrix(matG); // Simply duplicate matrix
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
    i_complex_matrix matLU_inv = initComplexMatrix(nSize, nSize);

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
i_complex_matrix pinv(const i_complex_matrix &matG)
{
    i_complex_matrix matGt = transpose(matG);
    i_complex_matrix matGtG_inv = inv(matMul(matGt, matG));
    return matMul(matGtG_inv, matGt);
}

// Moore-Penrose pseudoinversion (same as pinv(G) in MATLAB) [*1]
i_complex_matrix pinv2(const i_complex_matrix &matG, const i_float_t tolerance = 1.0e-9)
{
    bool useTranspose{false};
    const std::size_t nrows{matG.size()}, ncols{matG[0].size()};
    std::size_t nSize{ncols};

    i_complex_matrix matA, matGt;
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
        if (matA[i][i].real() > 0)
        {
            const i_float_t temp{std::abs(matA[i][i])};
            if (temp < tol)
            {
                tol = temp;
            }
        }
    }
    tol *= tolerance;

    i_complex_matrix matL = initComplexMatrix(nSize, nSize);
    std::size_t rankA{0};
    for (k = 0; k < nSize; ++k)
    {
        for (i = k; i < nSize; ++i)
        {
            matL[i][rankA] = matA[i][k];
            for (j = 0; j < rankA; ++j)
            {
                matL[i][rankA] -= matL[i][j] * std::conj(matL[k][j]);
            }
        }
        if (matL[k][rankA].real() > tol)
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
    i_complex_matrix matLt = transpose(matL);
    i_complex_matrix matM = inv(matMul(matLt, matL), false); // M = inv(L' * L)
    matA = matMul(matMul(matMul(matL, matM), matM), matLt);  // A = L * M * M * L'

    if (useTranspose)
    {
        return matMul(matGt, matA); // pinv(G) = G' * (L * M * M * L')
    }
    return matMul(matA, matGt); // pinv(G) = (L * M * M * L') * G'
}

// Calculate left division x = A \ b, using Moore-Penrose pinv, NOT same as MATLAB for a singular matrix
i_complex_matrix leftDiv(const i_complex_matrix &matA, const i_complex_matrix &matb)
{
    i_complex_matrix matx;
    if (matA.size() != matb.size())
    {
        std::cout << "Error when using leftDiv: row size not match.\n";
        return matx;
    }
    matx = matMul(pinv2(matA), matb); // x = A \ b = pinv(A) * b
    return matx;
}
