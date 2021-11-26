#include <bits/stdc++.h>
#include <unistd.h>
#ifdef __linux__
#define LIMPIAR "clear"
#endif // __linux__

#ifdef __MINGW32__
#define LIMPIAR "CLS"
#endif // __MINGW32__

/*Matrix with vectors*/
template <class t>
void print(std::vector<std::vector<t>> Mat)
{

    if (Mat.empty())
        std::cout << "La Matriz es de tamano 0\n";

    for (size_t i = 0; i < Mat.size(); i++)
    {
        for (size_t j = 0; j < Mat[i].size(); j++)
        {

            std::cout << Mat[i][j] << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

template <class t>
std::vector<std::vector<t>> resize(std::vector<std::vector<t>> Mat, int F, int C)
{
    std::vector<t> init(C, 0);
    std::vector<std::vector<t>> M(F, init);

    for (size_t i = 0; i < Mat.size(); i++)
    {
        for (size_t j = 0; j < Mat[i].size(); i++)
        {
            if (M[i][j] && Mat[i][j])
                M[i][j] = Mat[i][j];
        }
    }

    return M;
}

template <class t>
std::vector<std::vector<t>> ExtractMat(std::vector<std::vector<t>> Mat, int sz, int F, int C)
{

    std::vector<t> init(sz - 1, 0);
    std::vector<std::vector<t>> M(sz - 1, init);
    int k = 0, l = 0;

    for (int i = 0; i < sz; i++)
    {
        l = 0;
        for (int j = 0; j < sz; j++)
        {
            if (i != F && j != C)
            { //Si no estamos en la fila y la columna que se van a eliminar asignamos
                //el valor en esa posicion al valor k,l de la matriz resultado
                M[k][l] = Mat[i][j];
                l++; //iteramos las columnas de la matriz resultado
            }
        }
        if (i != F) //si no estamos en la fila que se va a eliminar iteramos k
            k++;    //iteramos las filas de la matriz resultado
    }

    return M;
}

template <class t>
t Det(std::vector<std::vector<t>> Mat, int sz)
{

    t detval;
    t dt;
    std::vector<t> v;

    detval = 0;

    if (sz == 1)
        return Mat[0][0];

    if (sz == 2)
    {
        detval = ((Mat[0][0] * Mat[1][1]) - (Mat[0][1] * Mat[1][0]));
        return detval;
    }
    else
    {
        std::vector<t> init(sz - 1, 0);
        std::vector<std::vector<t>> M(sz - 1, init);

        for (int i = 0; i < sz; i++)
        {
            M = ExtractMat(Mat, sz, 0, i);
            // res.print();

            dt = Det(M, sz - 1);

            //std::cout << "det -> " << dt << std::endl;

            if (i % 2 == 0)
            {
                v.push_back(Mat[0][i] * dt);
                //std::cout << Mat[0][i] * dt << std::endl;
                //std::cout << "+" << Mat[0][i] << "*" << dt << std::endl;
            }
            else
            {
                t men;
                men = -1;

                v.push_back((men * (Mat[0][i] * dt)));
                //std::cout << (men * (Mat[0][i] * dt)) << std::endl;
                //std::cout << "-(" << Mat[0][i] << "*" << dt << ")" << std::endl;
            }

            //std::cout << "detval = " << detval << std::endl;
        }

        for (int k = 0; k < v.size(); k++)
            detval = detval + v[k];
        return detval;
    }
}

template <class t>
t Determinante(std::vector<std::vector<t>> Mat)
{
    if (Mat.size() != Mat[0].size())
        throw std::invalid_argument("La matriz no es cuadrada");

    if (Mat.size() == 1 && Mat[0].size() == 1)
        return Mat[0][0];

    return Det(Mat, Mat.size());
}

template <class t>
std::vector<std::vector<t>> transp(std::vector<std::vector<t>> Mat)
{

    //std::cout << "DEBUG " << std::endl;

    std::vector<t> aux(Mat.size(), 0);
    std::vector<std::vector<t>> Trr(Mat[0].size(), aux);
    for (int i = 0; i < Mat[0].size(); i++)
    {
        for (int j = 0; j < Mat.size(); j++)
        {
            if (Mat[j][i])
                Trr[i][j] = Mat[j][i];
        }
    }

    return Trr;
}

template <class t>
std::vector<std::vector<t>> Adjunta(std::vector<std::vector<t>> Mat)
{

    //std::cout << "DEBUG " << std::endl;

    std::vector<t> aux(Mat.size(), 0);
    std::vector<std::vector<t>> M(Mat[0].size(), aux);

    std::vector<t> x(Mat.size() - 1, 0);
    std::vector<std::vector<t>> a(Mat[0].size() - 1, aux);

    for (int i = 0; i < Mat.size(); i++)
    {
        for (int j = 0; j < Mat[0].size(); j++)
        {
            a = ExtractMat(Mat, Mat.size(), i, j);

            M[i][j] = pow(-1, (i + 1) + (j + 1)) * Determinante(a);
        }
    }

    return M;
}
template <class t>
std::vector<std::vector<t>> Inversa(std::vector<std::vector<t>> Mat)
{

    //std::cout << "DEBUG " << std::endl;

    std::vector<t> aux(Mat.size(), 0);
    std::vector<std::vector<t>> M(Mat[0].size(), aux);

    std::vector<t> x(Mat.size() - 1, 0);
    std::vector<std::vector<t>> a(Mat[0].size() - 1, aux);

    a = Adjunta(Mat);

    M = transp(a);

    t d = Determinante(Mat);

    for (int i = 0; i < Mat.size(); i++)
    {
        for (int j = 0; j < Mat[0].size(); j++)
        {

            M[i][j] = M[i][j] / d;
        }
    }

    return M;
}

template <class t>
//Mat1 + Mat2
std::vector<std::vector<t>> Sum(std::vector<std::vector<t>> Mat1, std::vector<std::vector<t>> Mat2)
{
    if (Mat1.size() != Mat2.size() || Mat1[0].size() != Mat2[0].size())
        throw std::invalid_argument("Las Matrices no tienen el mismo tamano\n");

    std::vector<t> init(Mat1[0].size(), 5);
    std::vector<std::vector<t>> M(Mat1.size(), init);

    for (int i = 0; i < Mat1.size(); i++)
    {
        for (int j = 0; j < Mat1[0].size(); j++)
        {
            M[i][j] = Mat1[i][j] + Mat2[i][j];
        }
    }

    return M;
}
template <class t>
//Mat1-Mat2
std::vector<std::vector<t>> Res(std::vector<std::vector<t>> Mat1, std::vector<std::vector<t>> Mat2)
{
    if (Mat1.size() != Mat2.size() || Mat1[0].size() != Mat2[0].size())
        throw std::invalid_argument("Las Matrices no tienen el mismo tamano\n");

    std::vector<t> init(Mat1[0].size(), 5);
    std::vector<std::vector<t>> M(Mat1.size(), init);

    for (int i = 0; i < Mat1.size(); i++)
    {
        for (int j = 0; j < Mat1[0].size(); j++)
        {
            M[i][j] = Mat1[i][j] - Mat2[i][j];
        }
    }

    return M;
}

template <class t>
//Mat1-Mat2
std::vector<std::vector<t>> Mult(std::vector<std::vector<t>> Mat1, std::vector<std::vector<t>> Mat2)
{
    if (Mat1[0].size() != Mat2.size())
        throw std::invalid_argument("Las Matrices no se pueden multiplicar\n");

    std::vector<t> init(Mat2[0].size(), 0);
    std::vector<std::vector<t>> M(Mat1.size(), init);

    t sum; //si se planea usar objetos estos tienen que tener un 0 y sobrecargar el operador = para poder asignarlo

    sum = 0;

    for (size_t i = 0; i < Mat1.size(); i++)
    {
        sum = 0;
        for (size_t j = 0; j < Mat2[0].size(); j++)
        {
            for (size_t k = 0; k < Mat1[0].size(); k++)
            {

                M[i][j] += (Mat1[i][k] * Mat2[k][j]);
            }
        }
    }

    return M;
}

template <class t>
//Mat1-Mat2
std::vector<std::vector<t>> Multescalar(std::vector<std::vector<t>> Mat1, t alfa)
{

    std::vector<t> init(Mat1[0].size(), 0);
    std::vector<std::vector<t>> M(Mat1.size(), init);

    for (size_t i = 0; i < Mat1.size(); i++)
    {
        for (size_t j = 0; j < Mat1[0].size(); j++)
        {
            M[i][j] = alfa * Mat1[i][j];
        }
    }

    return M;
}

template <class t>
std::vector<std::vector<t>> extrac_partition(std::vector<std::vector<t>> Mat, size_t fila_init, size_t column_init, size_t fila_End, size_t column_end)
{
    int fi = ((fila_End - fila_init) + 1);
    int col = ((column_end - column_init) + 1);
    std::vector<t> init(col, 0);
    std::vector<std::vector<t>> M(fi, init);

    //Res.print();

    //std::cout << "filas: " << fi << "col: " << col << std::endl;

    int tr = 0;
    for (int i = fila_init; i <= fila_End; i++)
    {
        int k = 0;
        for (int j = column_init; j <= column_end; j++)
        {
            M[tr][k] = Mat[i][j];
            //if(i >= fila_init && j >= column_init)
            //std::cout << "->" << Mat[i][j] << '\t';
            k++;
        }
        //std::cout << std::endl;
        tr++;
    }

    return M;
}

template <class Ty>
std::vector<std::vector<Ty>> Inverse_partition(std::vector<std::vector<Ty>> Coef_Matriz, std::vector<std::vector<Ty>> vector_indp, int fila_init, int column_init)
{
    if (fila_init == Coef_Matriz.size() - 1)
    {
        std::cout << "No se pude obtener la inversa de la matriz , no hay mas particiones posibles!" << std::endl;
        return Coef_Matriz;
    }
    if (Coef_Matriz.size() != Coef_Matriz[0].size())
    {
        std::cout << "La matriz no es cuadrada!" << std::endl;
        return Coef_Matriz;
    }
    float det;

    std::vector<std::vector<Ty>> A11, A12, A21, A22, B1, B2;
    std::vector<std::vector<Ty>> A22inv, D, C, E, F, x1, x2, Cinv;

    std::vector<std::vector<Ty>> aux1, aux2, aux3, aux4, aux5, aux6, aux7, aux9, aux10, aux11;
    A22 = extrac_partition(Coef_Matriz, fila_init, column_init, Coef_Matriz.size() - 1, Coef_Matriz[0].size() - 1);
    A11 = extrac_partition(Coef_Matriz, 0, 0, fila_init - 1, column_init - 1);
    A12 = extrac_partition(Coef_Matriz, 0, column_init, fila_init - 1, Coef_Matriz[0].size() - 1);
    A21 = extrac_partition(Coef_Matriz, fila_init, 0, Coef_Matriz.size() - 1, column_init - 1);

    B1 = extrac_partition(vector_indp, 0, 0, fila_init - 1, 0);
    B2 = extrac_partition(vector_indp, fila_init, 0, vector_indp.size() - 1, 0);

    std::cout << "PARTICIONES: " << std::endl;

    std::cout << "A11" << std::endl;
    print<Ty>(A11);

    std::cout << "A12" << std::endl;
    print<Ty>(A12);

    std::cout << "A21" << std::endl;
    print<Ty>(A21);
    det = Determinante(A22);
    std::cout << "A22\n   Determinante  =  " << det << std::endl;
    print<Ty>(A22);

    std::cout << "B1 " << std::endl;
    print<Ty>(B1);

    std::cout << "B2" << std::endl;
    print<Ty>(B2);

    if (det != 0)
        A22inv = Inversa(A22);
    else
    {
        Inverse_partition(Coef_Matriz, vector_indp, fila_init + 1, column_init + 1);
    }
    std::cout << "A22^-1" << std::endl;
    print<Ty>(A22inv);

    try
    {
        D = Mult<Ty>(A12, A22inv);

        aux1 = Mult(D, A21);

        std::cout << "D * A21" << std::endl;
        print<Ty>(aux1);

        C = Res<Ty>(A11, aux1);

        Cinv = Inversa<Ty>(C);

        aux2 = Mult<Ty>(D, B2);
        aux3 = Res<Ty>(B1, aux2);
        x1 = Mult<Ty>(Cinv, aux3);

        E = Mult<Ty>(A22inv, A21);

        aux4 = Mult(E, Cinv);
        aux5 = Mult(aux4, D);
        F = Sum<Ty>(A22inv, aux5);

        aux6 = Mult<Ty>(F, B2);
        aux7 = Mult<Ty>(E, Cinv);
        aux9 = Mult<Ty>(aux7, B1);
        aux10 = Multescalar<Ty>(aux9, -1);

        x2 = Sum<Ty>(aux10, aux6);
    }
    catch (const std::exception &e)
    {
        std::system(LIMPIAR);
        Inverse_partition(Coef_Matriz, vector_indp, fila_init + 1, column_init + 1);
    }

    std::cout << "D" << std::endl;
    print<Ty>(D);

    std::cout << "C" << std::endl;
    print<Ty>(C);

    std::cout << "C^-1" << std::endl;
    print<Ty>(Cinv);

    std::cout << "E" << std::endl;
    print<Ty>(E);

    std::cout << "F" << std::endl;
    print<Ty>(F);

    std::cout << "X1" << std::endl;
    print<Ty>(x1);

    std::cout << "X2" << std::endl;
    print<Ty>(x2);

    return Coef_Matriz;
}

template <class Ty>
std::vector<std::vector<Ty>> Gauss_partition(std::vector<std::vector<Ty>> Coef_Matriz, std::vector<std::vector<Ty>> vector_indp, int fila_init, int column_init)
{
    if (fila_init == Coef_Matriz.size() - 1)
    {
        std::cout << "No se pude obtener la inversa de la matriz , no hay mas particiones posibles!" << std::endl;
        return Coef_Matriz;
    }
    if (Coef_Matriz.size() != Coef_Matriz[0].size())
    {
        std::cout << "La matriz no es cuadrada!" << std::endl;
        return Coef_Matriz;
    }
    float det, dt;

    std::vector<std::vector<Ty>> A11, A12, A21, A22, B1, B2;
    std::vector<std::vector<Ty>> A22pinv, A11inv, A12p, B1p, A22p, B2p, B2pp, B1pp;

    std::vector<std::vector<Ty>> aux1, aux2, aux3, aux4, aux5, aux6, aux7, aux9, aux10, aux11;
    A22 = extrac_partition(Coef_Matriz, fila_init, column_init, Coef_Matriz.size() - 1, Coef_Matriz[0].size() - 1);
    A11 = extrac_partition(Coef_Matriz, 0, 0, fila_init - 1, column_init - 1);
    A12 = extrac_partition(Coef_Matriz, 0, column_init, fila_init - 1, Coef_Matriz[0].size() - 1);
    A21 = extrac_partition(Coef_Matriz, fila_init, 0, Coef_Matriz.size() - 1, column_init - 1);

    B1 = extrac_partition(vector_indp, 0, 0, fila_init - 1, 0);
    B2 = extrac_partition(vector_indp, fila_init, 0, vector_indp.size() - 1, 0);

    std::cout << "PARTICIONES: " << std::endl;

    std::cout << "A11" << std::endl;
    dt = Determinante(A11);
    std::cout << "\nDeterminante  =  " << dt << std::endl;
    print<Ty>(A11);

    std::cout << "A12" << std::endl;
    print<Ty>(A12);

    std::cout << "A21" << std::endl;
    print<Ty>(A21);

    std::cout << "A22" << det << std::endl;
    print<Ty>(A22);

    std::cout << "B1 " << std::endl;
    print<Ty>(B1);

    std::cout << "B2" << std::endl;
    print<Ty>(B2);

    if (det != 0)
        A11inv = Inversa(A11);
    else
    {
        Inverse_partition(Coef_Matriz, vector_indp, fila_init + 1, column_init + 1);
    }

    std::cout << "A11^-1" << std::endl;
    print<Ty>(A11inv);

    try
    {
        A12p = Mult<Ty>(A11inv, A12);
        B1p = Mult<Ty>(A11inv, B1);
        aux1 = Mult<Ty>(A21, A12p);
        A22p = Res<Ty>(A22, aux1);

        aux2 = Mult<Ty>(A21, B1p);
        B2p = Res<Ty>(B2, aux2);

        if (Determinante(A22p) != 0)
            A22pinv = Inversa(A22p);
        else
        {
            Inverse_partition(Coef_Matriz, vector_indp, fila_init + 1, column_init + 1);
        }

        B2pp = Mult<Ty>(A22pinv, B2p);

        aux3 = Mult<Ty>(A12p, B2pp);
        B1pp = Res<Ty>(B1p, aux3);
    }
    catch (const std::exception &e)
    {
        std::system(LIMPIAR);
        Inverse_partition(Coef_Matriz, vector_indp, fila_init + 1, column_init + 1);
    }

    std::cout << "A12'" << std::endl;
    print<Ty>(A12p);

    std::cout << "B1'" << std::endl;
    print<Ty>(B1p);

    std::cout << "A11'" << std::endl;
    print<Ty>(A22p);

    std::cout << "B2'" << std::endl;
    print<Ty>(B2p);

    std::cout << "(A22')^-1" << std::endl;
    print<Ty>(A22pinv);

    std::cout << "B2''" << std::endl;
    print<Ty>(B2pp);

    std::cout << "B1''" << std::endl;
    print<Ty>(B1pp);

    return Coef_Matriz;
}

template <class Ty>
bool is_diagonaldominat(std::vector<std::vector<Ty>> Coef_Matriz)
{
    bool flag = true;
    for (size_t i = 0; i < Coef_Matriz.size(); i++)
    {
        Ty sum = 0;
        for (size_t j = 0; j < Coef_Matriz[0].size(); j++)
        {
            if (i != j)
                sum += Coef_Matriz[i][j];
        }

        if (std::abs(Coef_Matriz[i][i]) < std::abs(sum))
        {
            flag = false;
            break;
        }
    }
    return flag;
}

template <class Ty>
std::vector<std::vector<Ty>> Jacobi(std::vector<std::vector<Ty>> Coef_Matriz, std::vector<std::vector<Ty>> vector_indp, float error)
{
    if (!is_diagonaldominat(Coef_Matriz))
    {
        std::cout << "La matriz no es dominante en sentido diagonal " << std::endl;
        return Coef_Matriz;
    }

    std::vector<Ty> init(vector_indp[0].size(), 0);
    std::vector<std::vector<Ty>> X0(vector_indp.size(), init);

    std::vector<Ty> ini(vector_indp[0].size(), 0);
    std::vector<std::vector<Ty>> X1(vector_indp.size(), init);

    for (int i = 0; i < vector_indp.size(); i++)
        std::cout << "X" << i + 1 << "\t";
    std::cout << std::endl;

    print<Ty>(transp(X0));

    while (1)
    {
        for (size_t i = 0; i < vector_indp.size(); i++)
        {
            X1[i][0] = (1 / Coef_Matriz[i][i]);

            Ty sum = 0;

            for (size_t j = 0; j < vector_indp.size(); j++)
            {
                if (j != i)
                    sum += Coef_Matriz[i][j] * X0[j][0];
            }

            X1[i][0] = X1[i][0] * (vector_indp[i][0] - sum);
        }

        print<Ty>(transp(X1));

        bool flag = true;

        for (size_t i = 0; i < vector_indp.size(); i++)
        {
            if (std::abs(X1[i][0] - X0[i][0]) > error)
                flag = false;
        }
        X0 = X1;
        if (flag)
            break;
    }

    return vector_indp;
}

int main(int argc, char const *argv[])
{
    std::vector<float> init(4, 5);
    std::vector<std::vector<float>> M(4, init), X;

    for (size_t i = 0; i < 4; i++)
    {
        for (size_t j = 0; j < 4; j++)
        {
            //std::cout << "Ingrese " << i << j << ": ";
            std::cin >> M[i][j];
        }
    }

    std::vector<float> a(4, 0);
    std::vector<std::vector<float>> O(1, init), t;

    std::cout << "INdp" << std::endl;

    for (size_t i = 0; i < 1; i++)
    {
        for (size_t j = 0; j < 4; j++)
        {
            //std::cout << "Ingrese " << i << j << ": ";
            std::cin >> O[i][j];
        }
    }

    t = transp(O);
    X = Inversa(M);
    print<float>(M);
    print<float>(t);
    print<float>(X);

    Jacobi(M, t, 0.00001);

    return 0;
}
