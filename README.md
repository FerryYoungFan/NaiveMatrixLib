# NaiveMatrixLib 帆帆的简易矩阵计算库
A simple C++ stdlib-based complex &amp; real matrix library, with matrix inversion, left division (A\b) and determinant calculation.<br />
这是一个使用 C++ 标准库实现的简易复数及实数矩阵库，包含求逆矩阵、反斜杠除法 (A\b) 及行列式计算。

## Features 特点
* Designed for users who don't want to use large linear algebra libs.
* Only used C++ standard library, easy to learn and modify (Each file less than 600 lines).
* Header files only, separated complex and real matrix library.
* No recursive algorithm (using LU and Cholesky decomposition). Reliable for 1000 x 1000 and larger matrices.
- 如果你不想使用大型线性代数库来计算这些，那你来对地方了。（什么，你只是想交作业？）
- 仅使用C++标准库，无论是学习思维还是修改都很简单（每份代码都少于600行）。
- 仅使用头文件即可，复数和实数矩阵库是分开的。
- 没有递归运算（基于 LU 和 Cholesky 分解）。可对1000x1000及更大的矩阵使用。


## Available Functions 可用函数
* <b>rank</b>:    Matrix rank (Cholesky decomposition)
* <b>det</b>:     Matrix determinant calculation
* <b>inv</b>:     LU decomposition-based matrix inversion
* <b>pinv</b>:    pinv(G) = inv(G' * G) * G' (<b>WARNING</b>: full-rank matrix only!)
* <b>pinv2</b>:   Moore-Penrose pseudoinversion (same as pinv(G) in MATLAB)
* <b>leftDiv</b>: x = A \ b, using Moore-Penrose pinv, NOT same as MATLAB for a singular matrix
- <b>rank</b>:    计算矩阵秩 (Cholesky 分解)
- <b>det</b>:     计算矩阵行列式
- <b>inv</b>:     求逆矩阵，基于 LU 分解
- <b>pinv</b>:    经典伪逆，pinv(G) = inv(G' * G) * G' (<b>警告</b>：只能用于满秩矩阵！)
- <b>pinv2</b>:   Moore-Penrose 伪逆 (与 MATLAB 中的 pinv(G) 相同)
- <b>leftDiv</b>: 反斜杠除法 x = A \ b, 使用 Moore-Penrose 伪逆, 对奇异矩阵的处理与 MATLAB 不同


## Examples (Complex Matrices Only) 用例（仅列举复数矩阵）
```cpp
#include "matBasic_complex.hpp"
int main(int argc, char **argv)
{
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

    i_complex_matrix matb = {
        {{1.0, 0.0}},
        {{1.0, 2.0}},
        {{0.0, 3.0}}};
    showMatrix(matb, "b");
    showMatrix(leftDiv(matA, matb), "A \\ b");
    showMatrix(matMul(inv(matA), matb), "inv(A) * b");
    showMatrix(matMul(pinv(matA), matb), "pinv(A) * b");
    showMatrix(matMul(pinv2(matA), matb), "pinv2(A) * b (== A \\ b)");
    std::cin.get();
    return 0;
}
```
### Output 输出
```
A : 3 x 3 Complex Matrix:
    row[1]: 1,  0+1i,  4;
    row[2]: 0+5i,  1,  4;
    row[3]: 8,  1+1i,  0;

rank(A) = 3
det(A) = (-56,48)
inv(A) : 3 x 3 Complex Matrix:
    row[1]: 0.00588235+0.0764706i,  -0.00588235-0.0764706i,  0.0764706-0.00588235i;
    row[2]: -0.329412-0.282353i,  0.329412+0.282353i,  0.217647-0.170588i;
    row[3]: 0.177941+0.0632353i,  0.0720588-0.0632353i,  -0.0617647-0.0529412i;

pinv(A) : 3 x 3 Complex Matrix:
    row[1]: 0.00588235+0.0764706i,  -0.00588235-0.0764706i,  0.0764706-0.00588235i;
    row[2]: -0.329412-0.282353i,  0.329412+0.282353i,  0.217647-0.170588i;
    row[3]: 0.177941+0.0632353i,  0.0720588-0.0632353i,  -0.0617647-0.0529412i;

pinv2(A) : 3 x 3 Complex Matrix:
    row[1]: 0.00588235+0.0764706i,  -0.00588235-0.0764706i,  0.0764706-0.00588235i;
    row[2]: -0.329412-0.282353i,  0.329412+0.282353i,  0.217647-0.170588i;
    row[3]: 0.177941+0.0632353i,  0.0720588-0.0632353i,  -0.0617647-0.0529412i;



b : 3 x 1 Complex Matrix:
    row[1]: 1;
    row[2]: 1+2i;
    row[3]: 0+3i;

A \ b : 3 x 1 Complex Matrix:
    row[1]: 0.170588+0.217647i;
    row[2]: -0.0529412+1.31176i;
    row[3]: 0.535294-0.0411765i;

inv(A) * b : 3 x 1 Complex Matrix:
    row[1]: 0.170588+0.217647i;
    row[2]: -0.0529412+1.31176i;
    row[3]: 0.535294-0.0411765i;

pinv(A) * b : 3 x 1 Complex Matrix:
    row[1]: 0.170588+0.217647i;
    row[2]: -0.0529412+1.31176i;
    row[3]: 0.535294-0.0411765i;

pinv2(A) * b (== A \ b) : 3 x 1 Complex Matrix:
    row[1]: 0.170588+0.217647i;
    row[2]: -0.0529412+1.31176i;
    row[3]: 0.535294-0.0411765i;
```

## Reference 参考资料
* Pierre Courrieu, Fast Computation of Moore-Penrose Inverse Matrices, https://arxiv.org/abs/0804.4809
* Permute Sign Calculation, page5 https://www.math.rutgers.edu/docman-lister/math-main/academics/course-materials/250/assignments/1493-250c-lab3-sakai-pdf/file
* LU Decomposition C++ Implementation, https://blog.csdn.net/xx_123_1_rj/article/details/39553809
* LU Decomposition, https://www.math.ucdavis.edu/~linear/old/notes11.pdf
<br /><br /><br />
<p align="center">*** Project by Fanseline in Ericsson ***</p>
