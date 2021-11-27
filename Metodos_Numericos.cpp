#include <bits/stdc++.h>
#include <unistd.h>
#ifdef __linux__
#define LIMPIAR "clear"
#endif // __linux__

#ifdef __MINGW32__
#define LIMPIAR "CLS"
#endif // __MINGW32__
/*


*/

class EquationS
{

private:
    int equation;
    int method;

    inline float evaluar(float value, int func, bool derivate); // Evalua una funcion de acuerdo a un valor
    void bisection_method(std::pair<float, float> &interval, float error, float &root, bool print, int k, bool fil, std::ofstream &FILE);

    void falspos_method(std::pair<float, float> &interval, float error, float &root, bool print, int k, bool fil, std::ofstream &FILE);

    bool Newton_method(float n, float error, float &root, bool print, int k, bool fil, std::ofstream &FILE);

    void Secante_method(float n, float n_1, float error, float &root, bool print, int k, bool fil, std::ofstream &FILE);

public:
    EquationS(int method, int eq);

    float BiseccionM(std::pair<float, float> interval, float error, bool print, bool file);
    float FalsaPosicionM(std::pair<float, float> interval, float error, bool print, bool file);
    float NewtonM(float n, float error, bool print, bool file);

    float SecanteM(float n, float n_1, float error, bool print, bool file);

    ~EquationS();
};

EquationS::EquationS(int method, int eq)
{
    this->method = method;
    this->equation = eq;
}

inline float EquationS::evaluar(float value, int func, bool derivate)
{
    if (func == 1)
    {
        if (derivate) //retorna la evaluacion de la derivada
            return (-(sin(value)) - (cos(value)));
        else
            return (-(sin(value)) + cos(value));
    }
    else if (func == 2)
    {
        if (derivate) //retorna la evaluacion de la derivada
            return (exp(value) * (sin(value) + cos(value)));
        else
            return (exp(value) * sin(value));
    }
    else if (func == 3)
    {
        if (derivate) //retorna la evaluacion de la derivada
            return (pow(cos(value), 2) * ((2 * value) * (cos(value) - 3 * (pow(value, 2) - 5) * sin(value))));
        else
            return ((pow(value, 2) - 5) * pow(cos(value), 3));
    }
    else if (func == 4)
    {
        if (derivate) //retorna la evaluacion de la derivada
            return (sin(value) * (sin(value) + 2 * (value + 3) * cos(value)));
        else
            return ((value + 3) * pow(sin(value), 2));
    }

    return 0.0;
}

void EquationS::bisection_method(std::pair<float, float> &interval, float error, float &root, bool print, int k, bool fil, std::ofstream &FILE)
{

    float x = ((interval.first + interval.second) / 2.0);
    float fb = evaluar(interval.second, this->equation, false); // f(b)
    float fa = evaluar(interval.first, this->equation, false);  // f(a)
    float fx = evaluar(x, this->equation, false);

    if (std::abs(fa) < error)
    {
        root = interval.first;
        std::cout << "ITERACIONES: " << k << std::endl;
        return;
    }

    if (std::abs(fb) < error)
    {
        root = interval.second;
        std::cout << "ITERACIONES: " << k << std::endl;
        return;
    }

    if (print)
    {
        std::cout << "| " << interval.first << " | " << interval.second << " | ";
        std::cout << fa << " | " << fb << " | ";
        std::cout << x << " | " << fx << " | ";
        if (std::abs(fx) < error)
            std::cout << "TRUE";
        else
            std::cout << "FALSE";
        std::cout << " | ";
        if (fx * fa < 0)
            std::cout << "TRUE";
        else
            std::cout << "FALSE";
        std::cout << " | " << std::endl;
    }

    if (fil)
    {
        FILE << interval.first << ",";
        FILE << interval.second << ",";
        FILE << fa << "," << fb << ",";
        FILE << x << "," << fx << ",";

        if (std::abs(fx) < error)
            FILE << "TRUE";
        else
            FILE << "FALSE";

        FILE << ",";

        if (fx * fa < 0)
            FILE << "TRUE";
        else
            FILE << "FALSE";

        FILE << "\n";
    }

    if (std::abs(fx) < error)
    {
        root = x;
        std::cout << "ITERACIONES: " << k << std::endl;
        return;
    }

    std::pair<float, float> ninter;

    if (fx * fa < 0)
    {
        ninter = std::make_pair(interval.first, x);
    }
    else
    {
        ninter = std::make_pair(x, interval.second);
    }

    this->bisection_method(ninter, error, root, print, k + 1, fil, FILE);
    return;
}

float EquationS::BiseccionM(std::pair<float, float> interval, float error, bool print, bool file)
{
    float fb = this->evaluar(interval.second, this->equation, false); // f(b)
    float fa = this->evaluar(interval.first, this->equation, false);  // f(a)
    float dfb = this->evaluar(interval.second, this->equation, true); // f'(b)
    float dfa = this->evaluar(interval.first, this->equation, true);  // f'(a)
    float root = 0.0;

    std::ofstream File;
    File.open("Biseccion.csv");
    if (File.fail())
        perror("Error al crear el archivo! :(");

    std::cout << "f(" << interval.first << ") = " << fa << ", f(" << interval.second << ") = " << fb << std::endl;
    std::cout << "f'(" << interval.first << ") = " << dfa << ", f'(" << interval.second << ") = " << dfb << std::endl;

    if (!(fa * fb < 0 && dfa * dfb > 0))
    {
        std::cout << "NO HAY RAIZ AISLADA :( " << std::endl;
        return 0;
    }
    //throw std::invalid_argument(":("); //No hay raiz aislada,

    if (print)
        std::cout << "|  a   |  b  |   f(a)   |   f(b)   | x = (a+b)/2 |  f(x)  |   |f(x)| < " << error << " |   f(x)*f(a) < 0   | " << std::endl;

    if (file)
    {
        File << "a,b,f(a),f(b),x = (a+b)/2, f(x) ,|f(x)| < ";
        File << error;
        File << ",f(x)*f(a) < 0   \n";
    }

    this->bisection_method(interval, error, root, print, 0, file, File);

    std::cout << "LA RAIZ ES : " << root << " ( aprox. )" << std::endl;

    File.close();

    return root;
}

void EquationS::falspos_method(std::pair<float, float> &interval, float error, float &root, bool print, int k, bool fil, std::ofstream &FILE)
{

    float fb = evaluar(interval.second, this->equation, false); // f(b)
    float fa = evaluar(interval.first, this->equation, false);  // f(a)
    float x = (((interval.second * fa) - (interval.first * fb)) / (fa - fb));
    float fx = evaluar(x, this->equation, false);

    if (std::abs(fa) < error)
    {
        root = interval.first;
        std::cout << "ITERACIONES: " << k << std::endl;
        return;
    }

    if (std::abs(fb) < error)
    {
        root = interval.second;
        std::cout << "ITERACIONES: " << k << std::endl;
        return;
    }

    if (print)
    {
        std::cout << "| " << interval.first << " | " << interval.second << " | ";
        std::cout << fa << " | " << fb << " | ";
        std::cout << x << " | " << fx << " | ";
        if (std::abs(fx) < error)
            std::cout << "TRUE";
        else
            std::cout << "FALSE";
        std::cout << " | ";
        if (fx * fa < 0)
            std::cout << "TRUE";
        else
            std::cout << "FALSE";
        std::cout << " | " << std::endl;
    }

    if (fil)
    {
        FILE << interval.first << ",";
        FILE << interval.second << ",";
        FILE << fa << "," << fb << ",";
        FILE << x << "," << fx << ",";

        if (std::abs(fx) < error)
            FILE << "TRUE";
        else
            FILE << "FALSE";

        FILE << ",";

        if (fx * fa < 0)
            FILE << "TRUE";
        else
            FILE << "FALSE";

        FILE << "\n";
    }

    if (std::abs(fx) < error)
    {
        root = x;
        std::cout << "ITERACIONES: " << k << std::endl;
        return;
    }

    std::pair<float, float> ninter;

    if (fx * fa < 0)
    {
        ninter = std::make_pair(interval.first, x);
    }
    else
    {
        ninter = std::make_pair(x, interval.second);
    }

    this->bisection_method(ninter, error, root, print, k + 1, fil, FILE);
    return;
}

float EquationS::FalsaPosicionM(std::pair<float, float> interval, float error, bool print, bool file)
{
    float fb = this->evaluar(interval.second, this->equation, false); // f(b)
    float fa = this->evaluar(interval.first, this->equation, false);  // f(a)
    float dfb = this->evaluar(interval.second, this->equation, true); // f'(b)
    float dfa = this->evaluar(interval.first, this->equation, true);  // f'(a)
    float root = 0.0;

    std::ofstream File;
    File.open("FalsaPosicion.csv");
    if (File.fail())
        perror("Error al crear el archivo! :(");

    std::cout << "f(" << interval.first << ") = " << fa << ", f(" << interval.second << ") = " << fb << std::endl;
    std::cout << "f'(" << interval.first << ") = " << dfa << ", f'(" << interval.second << ") = " << dfb << std::endl;

    if (!(fa * fb < 0 && dfa * dfb > 0))
        return 0;
    //throw std::invalid_argument(":("); //No hay raiz aislada,

    if (print)
        std::cout << "|  a   |  b  |   f(a)   |   f(b)   | x = bf(a)-af(b)/ f(a)-f(b) |  f(x)  |   |f(x)| < " << error << " |   f(x)*f(a) < 0   | " << std::endl;

    if (file)
    {
        File << "a,b,f(a),f(b),x = bf(a)-af(b)/ f(a)-f(b) , f(x) ,|f(x)| < ";
        File << error;
        File << "f(x)*f(a) < 0   \n";
    }

    this->bisection_method(interval, error, root, print, 0, file, File);

    std::cout << "LA RAIZ ES : " << root << " ( aprox. )" << std::endl;

    return root;
}

bool EquationS::Newton_method(float n, float error, float &root, bool print, int k, bool fil, std::ofstream &FILE)
{
    float fn = this->evaluar(n, this->equation, false); // f(b)
    float dfn = this->evaluar(n, this->equation, true); // f(a)

    if (dfn == 0 || k > 1000)
        return false;

    float np1 = n - (fn / dfn);
    float fnp1 = this->evaluar(np1, this->equation, false);

    if (std::abs(fn) < error)
    {
        root = n;
        std::cout << "ITERACIONES: " << k << std::endl;
        return true;
    }

    if (print)
    {
        std::cout << "| " << n << " | ";
        std::cout << fn << " | " << dfn << " | ";
        std::cout << np1 << " | " << fnp1 << " | ";
        if (std::abs(fnp1) < error)
            std::cout << "TRUE";
        else
            std::cout << "FALSE";

        std::cout << " | " << std::endl;
    }

    if (fil)
    {

        FILE << n << ",";
        FILE << fn << "," << dfn << ",";
        FILE << np1 << "," << fnp1 << ",";

        if (std::abs(fnp1) < error)
            FILE << "TRUE";
        else
            FILE << "FALSE";

        FILE << "\n";
    }

    if (std::abs(fnp1) < error)
    {
        root = np1;
        std::cout << "ITERACIONES: " << k << std::endl;
        return true;
    }

    return this->Newton_method(np1, error, root, print, k + 1, fil, FILE);
}

float EquationS::NewtonM(float n, float error, bool print, bool file)
{
    float fn = this->evaluar(n, this->equation, false); // f(b)
    float dfn = this->evaluar(n, this->equation, true); // f(a)
    float root = 0.0;

    std::ofstream File;
    File.open("Metodo_De_Newton.csv");
    if (File.fail())
        perror("Error al crear el archivo! :(");

    if (dfn == 0)
    {
        std::cout << "Metodo indeterminado f'(x) = 0 !" << std::endl;
        return 0;
    }

    std::cout << "f(" << n << ") = " << fn << ", f'(" << n << ") = " << dfn << std::endl;

    if (print)
        std::cout << "| Xn |  f(Xn)   |   f'(Xn)   | Xn+1 = Xn - (f(Xn)/f'(Xn)) |  f(Xn+1)  |   |f(Xn+1)| < " << error << std::endl;

    if (file)
    {
        File << "Xn,f(Xn),f'(Xn),Xn+1 = Xn - (f(Xn)/f'(Xn)) , f(Xn+1) ,|f(Xn+1)| < ";
        File << error;
        File << "\n";
    }

    if (this->Newton_method(n, error, root, print, 0, file, File))
        std::cout << "LA RAIZ ES : " << root << " ( aprox. )" << std::endl;
    else
        std::cout << "Metodo indeterminado f'(x) = 0 !" << std::endl;

    return root;
}

void EquationS::Secante_method(float n, float n_1, float error, float &root, bool print, int k, bool fil, std::ofstream &FILE)
{
    float fn = this->evaluar(n, this->equation, false);
    float fn_1 = this->evaluar(n_1, this->equation, false);
    float xnp1;
    if (fn - fn_1 != 0)
        xnp1 = n - (((fn) * (n - n_1)) / (fn - fn_1));
    float fnp1 = this->evaluar(xnp1, this->equation, false);

    if (std::abs(fn) < error)
    {
        root = n;
        std::cout << "ITERACIONES: " << k << std::endl;
        return;
    }

    if (std::abs(fn_1) < error)
    {
        root = n_1;
        std::cout << "ITERACIONES: " << k << std::endl;
        return;
    }

    if (print)
    {
        std::cout << "| " << n_1 << " | ";
        std::cout << n << " | " << fn << " | ";
        std::cout << fn_1 << " | " << xnp1 << " | ";
        if (std::abs(fnp1) < error)
            std::cout << "TRUE";
        else
            std::cout << "FALSE";

        std::cout << " | " << std::endl;
    }

    if (fil)
    {

        FILE << n_1 << ",";
        FILE << n << "," << fn << ",";
        FILE << fn_1 << "," << xnp1 << ",";

        if (std::abs(fnp1) < error)
            FILE << "TRUE";
        else
            FILE << "FALSE";

        FILE << "\n";
    }

    if (std::abs(fnp1) < error)
    {
        root = xnp1;
        std::cout << "ITERACIONES: " << k << std::endl;
        return;
    }

    this->Secante_method(xnp1, n, error, root, print, k + 1, fil, FILE);
}

float EquationS::SecanteM(float n, float n_1, float error, bool print, bool file)
{
    float fn = this->evaluar(n, this->equation, false);     // f(b)
    float fn_1 = this->evaluar(n_1, this->equation, false); // f(a)
    float root = 0.0;

    std::ofstream File;
    File.open("Metodo_De_La_Secante.csv");
    if (File.fail())
        perror("Error al crear el archivo! :(");

    if (print)
        std::cout << "| Xn |  Xn-1   |   f(Xn)   | f(Xn-1) |  Xn+1 = Xn - ((f(Xn))*(Xn-(Xn-1))/(f(Xn)-f(Xn-1)))  |   |f(Xn+1)| < " << error << std::endl;

    if (file)
    {
        File << "Xn,f(Xn),f'(Xn),Xn+1 = Xn - (f(Xn)/f'(Xn)) , f(Xn+1) ,|f(Xn+1)| <";
        File << error;
        File << "\n";
    }

    this->Secante_method(n, n_1, error, root, print, 0, file, File);

    std::cout << "LA RAIZ ES : " << root << " ( aprox. )" << std::endl;

    return root;
}

EquationS::~EquationS()
{
}

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
    std::vector<std::vector<Ty>> X1(vector_indp.size(), init), Comp;

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
    Comp = Mult<Ty>(Coef_Matriz, X1);
    std::cout << "Comprobacion" << std::endl;
    print<Ty>(Comp);
    return X1;
}

template <class Ty>
std::vector<std::vector<Ty>> GaussSeidel(std::vector<std::vector<Ty>> Coef_Matriz, std::vector<std::vector<Ty>> vector_indp, float error)
{
    if (!is_diagonaldominat(Coef_Matriz))
    {
        std::cout << "La matriz no es dominante en sentido diagonal " << std::endl;
        return Coef_Matriz;
    }

    std::vector<Ty> init(vector_indp[0].size(), 0);
    std::vector<std::vector<Ty>> X0(vector_indp.size(), init);

    std::vector<Ty> ini(vector_indp[0].size(), 0);
    std::vector<std::vector<Ty>> X1(vector_indp.size(), init), Comp;

    for (int i = 0; i < vector_indp.size(); i++)
        std::cout << "X" << i + 1 << "\t";
    std::cout << std::endl;

    print<Ty>(transp(X0));

    while (1)
    {
        for (size_t i = 0; i < vector_indp.size(); i++)
        {
            X1[i][0] = (1 / Coef_Matriz[i][i]);

            Ty sum = 0, sum2 = 0;

            for (size_t j = 0; j < i; j++)
            {
                if (j != i)
                    sum += Coef_Matriz[i][j] * X1[j][0];
            }
            for (size_t k = i + 1; k < vector_indp.size(); k++)
            {
                if (k != i)
                    sum2 += Coef_Matriz[i][k] * X0[k][0];
            }

            X1[i][0] = X1[i][0] * (vector_indp[i][0] - sum - sum2);
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
    Comp = Mult<Ty>(Coef_Matriz, X1);

    std::cout << "Comprobacion" << std::endl;
    print<Ty>(Comp);
    return X1;
}

template <class Ty>
std::vector<std::vector<Ty>> Sustitucion_adelante(std::vector<std::vector<Ty>> L, std::vector<std::vector<Ty>> vector_indp)
{
    std::vector<Ty> init(vector_indp[0].size(), 0);
    std::vector<std::vector<Ty>> C(vector_indp.size(), init);

    for (size_t i = 0; i < vector_indp.size(); i++)
    {
        for (size_t j = 0; j < vector_indp.size(); j++)
        {
            C[i][0] = vector_indp[i][0] / L[i][i];

            if (i - 1 >= 0 && i != j)
                C[i][0] += -((C[j][0] * L[i][j]) / L[i][i]);
        }
    }

    return C;
}

template <class Ty>
std::vector<std::vector<Ty>> Sustitucion_atras(std::vector<std::vector<Ty>> U, std::vector<std::vector<Ty>> vector_indp, std::vector<std::vector<Ty>> C)
{
    std::vector<Ty> init(vector_indp[0].size(), 0);
    std::vector<std::vector<Ty>> X(vector_indp.size(), init);

    for (size_t i = vector_indp.size() - 1; i >= 0; i++)
    {
        for (size_t j = vector_indp.size() - 1; j >= 0; j++)
        {
            X[i][0] = C[i][0] / U[i][i];

            if (i + 1 <= vector_indp.size() - 1 && i != j)
                C[i][0] += -((C[j][0] * U[i][j]) / U[i][i]);
        }
    }

    return X;
}
template <class Ty>
bool is_positive(std::vector<std::vector<Ty>> U)
{
    bool f = true;

    for (size_t i = 1; i <= U.size(); i++)
    {
        std::vector<Ty> init(i, 0);
        std::vector<std::vector<Ty>> X(i, init);

        for (size_t j = 0; j < i; j++)
        {
            for (size_t k = 0; k < i; k++)
            {
                X[j][k] = U[j][k];
            }
        }

        if (Determinante(U) > 0)
            continue;
        else
        {
            f = false;
            break;
        }
    }

    return f;
}

template <class Ty>
bool is_simetric(std::vector<std::vector<Ty>> U)
{
    std::vector<std::vector<Ty>> T;
    bool f = true;

    T = transp(U);

    //print<Ty>(U);
    //print<Ty>(T);

    for (size_t i = 0; i < U.size(); i++)
    {
        for (size_t j = 0; j < U.size(); j++)
        {
            if (T[i][j] == U[i][j])
            {
                f = true;
            }
            else
            {
                f = false;
                break;
            }
        }
    }

    return f;
}

template <class Ty>
std::vector<std::vector<Ty>> cholesky(std::vector<std::vector<Ty>> U, std::vector<std::vector<Ty>> vector_indp)
{
    // Decomposing a matrix into Lower Triangular
    std::vector<Ty> init(U[0].size(), 0);
    std::vector<std::vector<Ty>> lower(U.size(), init);
    std::vector<std::vector<Ty>> ltranpose;

    //if (!is_positive(U))
    //{
    //    std::cout << "La matriz no esta definida positivamente!" << std::endl;
    //    return U;
    //}

    //if (!is_simetric(U))
    //{
    //    std::cout << "La matriz no es simetrica!" << std::endl;
    //    return U;
    //}

    for (int i = 0; i < U.size(); i++)
    {
        for (int j = 0; j <= i; j++)
        {
            Ty sum = 0;

            if (j == i) // summation for diagonals
            {
                for (int k = 0; k < j; k++)
                    sum += pow(lower[j][k], 2);
                lower[j][j] = sqrt(U[j][j] - sum);
            }
            else
            {

                // Evaluating L(i, j) using L(j, j)
                for (int k = 0; k < j; k++)
                    sum += (lower[i][k] * lower[j][k]);
                lower[i][j] = (U[i][j] - sum) / lower[j][j];
            }
        }
    }

    ltranpose = transp(lower);
    std::cout << "L " << std::endl;
    print<Ty>(lower);

    std::cout << "U " << std::endl;
    print<Ty>(ltranpose);

    return U;
}

template <class Ty>
std::vector<std::vector<Ty>> Doolittle(std::vector<std::vector<Ty>> mat, std::vector<std::vector<Ty>> vector_indp)
{
    // Decomposing a matrix into Lower Triangular
    std::vector<Ty> init(mat[0].size(), 0);
    std::vector<std::vector<Ty>> lower(mat.size(), init);
    std::vector<Ty> in(mat[0].size(), 0);
    std::vector<std::vector<Ty>> upper(mat.size(), init);

    for (int i = 0; i < mat.size(); i++)
    {
        // Upper Triangular
        for (int k = i; k < mat.size(); k++)
        {
            // Summation of L(i, j) * U(j, k)
            Ty sum = 0;
            for (int j = 0; j < i; j++)
                sum += (lower[i][j] * upper[j][k]);

            // Evaluating U(i, k)
            upper[i][k] = mat[i][k] - sum;
        }

        // Lower Triangular
        for (int k = i; k < mat.size(); k++)
        {
            if (i == k)
                lower[i][i] = 1; // Diagonal as 1
            else
            {
                // Summation of L(k, j) * U(j, i)
                Ty sum = 0;
                for (int j = 0; j < i; j++)
                    sum += (lower[k][j] * upper[j][i]);

                // Evaluating L(k, i)
                lower[k][i] = (mat[k][i] - sum) / upper[i][i];
            }
        }
    }

    std::cout << "L " << std::endl;
    print<Ty>(lower);

    std::cout << "U " << std::endl;
    print<Ty>(upper);

    return mat;
}

class Menu
{
private:
    int op_main_menu; //Opcion del menu principal
    /*
        Menu Principal:
            1. Solucion de ecuacionnes 
            2. Sistemas de ecuaciones lineales
            3. Factorizacion LU
    */
    int op_sub_menu1; //Sub menu de Soluciones de ecuaciones
    /*
        Sub Menu Soluciones de ecuaciones:
            1.1 Método de Bisección.
            1.2 Método de falsa posición.
            1.3 Método de Newton.
            1.4 Método de la secante.
    */

    int function_submenu1; // Sub menu de ecuaciones para el submenu 1

    int op_submenu_2; //Sub menu sistemas de ecuaciones
    int op_metodos_Exactos;
    int op_metodos_iter;

    int op_submenu3;

    std::vector<float> readNumeric_input(int argsN); //lee una linea y la retorna como un string
    bool read_input(std::string find);               // lee un input y busca una subcadena(find) retorna true si existe una ocurrencia

public:
    Menu();
    template <class F>
    std::vector<std::vector<F>> read_matrix(size_t filas, size_t columnas, bool n); //lee una matriz de n * m y retorna la matriz

    int Main_menu();

    void Menu_Soluciones_de_ecuaciones();
    int Menu_ecuaciones();
    void Menu_Biseccion(int equa);
    void Menu_FalsaPos(int equa);
    void Menu_Newton(int equa);
    void Menu_Secante(int equa);

    void Menu_sistemas_ecuaciones();
    void Menu_metodos_exactos();
    void Menu_inversa_particionado();
    void Menu_Gausss_JordanP();
    void Menu_Metodo_intercambio();

    void Menu_metodos_iterativos();

    void Menu_jacobi();
    void Menu_gaussseidel();
    void Menu_relajacion();

    void Menu_factorizacionLU();
    void Menu_cholesky();
    void Menu_Doolittle();

    ~Menu();
};

Menu::Menu()
{
    this->op_main_menu = 0;
    this->op_sub_menu1 = 0;
    this->function_submenu1 = 0;
}

std::vector<float> Menu::readNumeric_input(int argsN)
{ //Lee una entrada numerica con un numero n de argumentos y retorna un vector con dichos elementos
    //*si el usuario introduce mas de los N argumentos, dichos seran descartados
    std::string s = "";
    char line[100];
    std::vector<float> v;

    //std::cin.ignore(1024, '\n');
    //std::cout << "HOla" << std::endl;
    //std::cin
    std::cin.getline(line, 100, '\n');

    //std::cout << line << std::endl;

    for (size_t i = 0; line[i] && v.size() <= argsN; i++)
    {

        //std::cout << line[i] << std::endl;
        if ((line[i] - '0' >= 0 && line[i] - '0' <= 9) || (line[i] == '.') || (line[i] == '-'))
        {
            s += line[i];
        }

        //std::cout << "s-> " << s << std::endl;

        if (line[i] == ' ' || line[i] == ',')
        {

            try
            {
                v.push_back(std::stof(s));
            }
            catch (const std::exception &e)
            {
                std::cerr << e.what() << '\n';
            }

            s = "";
        }
    }

    if (v.size() == 0 || s != "")
    {
        try
        {
            v.push_back(std::stof(s));
        }
        catch (const std::exception &e)
        {
        }
    }

    //for (auto c : v)
    //std::cout << c << ' ';

    std::cout << std::endl;

    return v;
}

bool Menu::read_input(std::string find)
{
    std::string s;

    std::getline(std::cin, s);

    std::transform(s.begin(), s.end(), s.begin(), tolower);

    size_t f = s.find(find);

    if (f != std::string::npos)
        return true;

    return false;
}

template <class F>
std::vector<std::vector<F>> Menu::read_matrix(size_t filas, size_t columnas, bool n)
{
    std::vector<std::string> v;
    std::vector<std::vector<float>> m;
    std::vector<F> init(columnas, 0);
    std::vector<std::vector<F>> M(filas, init);

    for (size_t i = 0; i < filas; i++)
    {
        std::string s;
        std::getline(std::cin, s);
        v.push_back(s);
    }

    std::string aux = "", d = "", nu = "";
    for (int i = 0; i < filas; i++)
    {
        int k = 0;
        std::vector<float> x;
        for (int j = 0; j < v[i].size() && k < columnas; j++)
        {

            if ((v[i][j] - '0' >= 0 && v[i][j] - '0' <= 9) || (v[i][j] == '.') || (v[i][j] == '-') || (v[i][j] == '/'))
            {
                if (v[i][j] != '/')
                    aux += v[i][j];
                else
                {
                    nu = aux;
                    aux += v[i][j];

                    int t = j + 1;
                    while (t < v[i].size() && (v[i][t] != ' ' || v[i][t] != ',' || d == ""))
                    {
                        d += v[i][t];
                        t++;
                    }
                }
            }

            if (v[i][j] == ' ' || v[i][j] == ',' || j == v[i].size() - 1)
            {
                try
                {
                    //std::cout << "Insert: " << aux << std::endl;
                    if (aux.find('/') == std::string::npos)
                        M[i][k] = std::stof(aux);
                    else
                    {
                        float den, nume;
                        try
                        {
                            den = std::stof(d);
                            nume = std::stof(nu);
                        }
                        catch (const std::exception &e)
                        {
                            throw std::invalid_argument("invalid input frac!");
                        }

                        M[i][k] = nume / den;

                        nu = "";
                        d = "";
                    }
                    aux = "";
                    k++;
                }
                catch (const std::exception &e)
                {
                    throw std::invalid_argument("invalid input!");
                }
            }
        }
    }

    //std::cout << "matrix :" << std::endl;

    //Mat.print();

    return M;
}

int Menu::Main_menu()
{ //Menu Principal
    std::system(LIMPIAR);
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << "|                   METODOS NUMERICOS                      |" << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << "| 1 . SOLUCIONES DE ECUACIONES                             |" << std::endl;
    std::cout << "|----------------------------------------------------------|" << std::endl;
    std::cout << "| 2 . SISTEMAS DE ECUACIONES LINEALES                      |" << std::endl;
    std::cout << "|----------------------------------------------------------|" << std::endl;
    std::cout << "| 3 . FACTORIZACION LU                                     |" << std::endl;
    std::cout << "|----------------------------------------------------------|" << std::endl;
    std::cout << "| 4 . SALIR                                                |" << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;
    std::vector<float> v;
    do
    {

        std::cout << "Ingrese una opcion del menu(numero): "; // << std::endl;

        v = this->readNumeric_input(10);

        //std::cout << (v.size() > 0) ? v[0] : -1;
        //std::cout << "\n";

        if (v.size() >= 1 && ((((int)v[0]) >= 1) && (((int)v[0]) <= 4)))
        {
            break;
        }
    } while (1);
    //sleep(100);
    this->op_main_menu = v[0];

    //std::cout << this->op_main_menu << std::endl;
    //sleep(10);
    std::system(LIMPIAR);

    switch (this->op_main_menu)
    {
    case 1:

        this->Menu_Soluciones_de_ecuaciones();
        break;
    case 2:
        this->Menu_sistemas_ecuaciones();
    case 3:
        this->Menu_factorizacionLU();
        break;
    default:
        break;
    }

    return v[0];
}

void Menu::Menu_Soluciones_de_ecuaciones()
{
    do
    {
        std::system(LIMPIAR);
        std::cout << "------------------------------------------------------------" << std::endl;
        std::cout << "|                SOLUCIONES DE ECUACIONES                  |" << std::endl;
        std::cout << "------------------------------------------------------------" << std::endl;
        std::cout << "| 1 . METODO DE BISECCION                                  |" << std::endl;
        std::cout << "|----------------------------------------------------------|" << std::endl;
        std::cout << "| 2 . METODO DE FALSA POSICION                             |" << std::endl;
        std::cout << "|----------------------------------------------------------|" << std::endl;
        std::cout << "| 3 . METODO DE NEWTON                                     |" << std::endl;
        std::cout << "|----------------------------------------------------------|" << std::endl;
        std::cout << "| 4 . METODO DE LA SECANTE                                 |" << std::endl;
        std::cout << "------------------------------------------------------------" << std::endl;

        std::vector<float> v;
        int equa;
        do
        {

            std::cout << "Ingrese una opcion del menu(numero): "; // << std::endl;

            v = this->readNumeric_input(1);

            if (v.size() >= 1 && ((((int)v[0]) >= 1) && (((int)v[0]) <= 4)))
            {
                break;
            }
        } while (1);

        do
        {
            this->op_sub_menu1 = ((int)v[0]);

            equa = this->Menu_ecuaciones();

            switch (((int)v[0])) //Metodo seleccionado
            {
            case 1:
                this->Menu_Biseccion(equa);
                break;

            case 2:
                this->Menu_FalsaPos(equa);
                break;

            case 3:
                this->Menu_Newton(equa);
                break;

            case 4:
                this->Menu_Secante(equa);
                break;

            default:
                break;
            }

            std::cout << "¿Deseas probar otra ecuacion?(S/N) : ";
            if (!this->read_input("s"))
                break;
        } while (1);

        std::cout << "¿Deseas probar otro metodo?(S/N) : ";
        if (!this->read_input("s"))
            break;

    } while (1);

    //std::cout << op_sub_menu1 << std::endl;
}

int Menu::Menu_ecuaciones()
{
    std::vector<float> v;
    std::system(LIMPIAR);
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << "|                       ECUACIONES                         |" << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << "| 1 . f(x) = -SEN(X) + COS(X)                              |" << std::endl;
    std::cout << "|----------------------------------------------------------|" << std::endl;
    std::cout << "| 2 . f(x) =  e^x*SEN(X)                                   |" << std::endl;
    std::cout << "|----------------------------------------------------------|" << std::endl;
    std::cout << "| 3 . f(x) =  (x^2 - 5)*COS^3(X)                           |" << std::endl;
    std::cout << "|----------------------------------------------------------|" << std::endl;
    std::cout << "| 4 . f(x) =  (x+3)*SEN^2(X)                               |" << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;

    do
    {

        std::cout << "Ingrese una opcion del menu(numero): "; // << std::endl;

        v = this->readNumeric_input(1);

        if (v.size() >= 1 && ((((int)v[0]) >= 1) && (((int)v[0]) <= 4)))
        {
            break;
        }
    } while (1);

    //std::cout << "this:: " << v[0] << std::endl;
    //sleep(10);

    return ((int)v[0]);
}

void Menu::Menu_Biseccion(int equa)
{
    EquationS metod(1, equa);

    std::vector<float> interval, err;

    do
    {
        std::system(LIMPIAR);
        if (equa == 1)
        {
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "|                    BISECCION                             |" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "| 1 . f(x) = -SEN(X) + COS(X)                              |" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
        }
        else if (equa == 2)
        {
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "|                    BISECCION                             |" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "| 2 . f(x) =  e^x*SEN(X)                                   |" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
        }
        else if (equa == 3)
        {
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "|                    BISECCION                             |" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "| 3 . f(x) =  (x^2 - 5)*COS^3(X)                           |" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
        }
        else if (equa == 4)
        {
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "|                    BISECCION                             |" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "| 4 . f(x) =  (x+3)*SEN^2(X)                               |" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
        }

        do
        {
            std::cout << "Ingrese el intervalo(any format): ";

            interval = this->readNumeric_input(2);

            if (interval.size() >= 2)
            {
                break;
            }

        } while (1);

        do
        {
            std::cout << "Ingrese el error de tolerancia: ";
            err = this->readNumeric_input(1);

            if (err.size() > 0)
                break;
        } while (1);

        metod.BiseccionM(std::make_pair(interval[0], interval[1]), err[0], true, true);

        std::cout << "¿Deseas probar otro intervalo?(S/N) : ";
        if (!this->read_input("s"))
            break;

    } while (1);
}

void Menu::Menu_FalsaPos(int equa)
{
    EquationS metod(2, equa);

    std::vector<float> interval, err;

    do
    {
        std::system(LIMPIAR);
        if (equa == 1)
        {
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "|                    FALSA POSICION                        |" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "| 1 . f(x) = -SEN(X) + COS(X)                              |" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
        }
        else if (equa == 2)
        {
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "|                    FALSA POSICION                        |" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "| 2 . f(x) =  e^x*SEN(X)                                   |" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
        }
        else if (equa == 3)
        {
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "|                    FALSA POSICION                        |" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "| 3 . f(x) =  (x^2 - 5)*COS^3(X)                           |" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
        }
        else if (equa == 4)
        {
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "|                    FALSA POSICION                        |" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "| 4 . f(x) =  (x+3)*SEN^2(X)                               |" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
        }
        do
        {
            std::cout << "Ingrese el intervalo(any format): ";

            interval = this->readNumeric_input(2);

            if (interval.size() >= 2)
            {
                break;
            }

        } while (1);

        do
        {
            std::cout << "Ingrese el error de tolerancia: ";
            err = this->readNumeric_input(1);

            if (err.size() > 0)
                break;
        } while (1);

        metod.FalsaPosicionM(std::make_pair(interval[0], interval[1]), err[0], true, true);

        std::cout << "¿Deseas probar otro intervalo?(S/N) : ";
        if (!this->read_input("s"))
            break;

    } while (1);
}

void Menu::Menu_Newton(int equa)
{
    EquationS metod(3, equa);

    std::vector<float> aprox, err;

    do
    {
        std::system(LIMPIAR);
        if (equa == 1)
        {
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "|                    Newton Raphson                        |" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "| 1 . f(x) = -SEN(X) + COS(X)                              |" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
        }
        else if (equa == 2)
        {
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "|                    Newton Raphson                        |" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "| 2 . f(x) =  e^x*SEN(X)                                   |" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
        }
        else if (equa == 3)
        {
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "|                    Newton Raphson                        |" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "| 3 . f(x) =  (x^2 - 5)*COS^3(X)                           |" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
        }
        else if (equa == 4)
        {
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "|                    Newton Raphson                        |" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "| 4 . f(x) =  (x+3)*SEN^2(X)                               |" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
        }
        do
        {
            std::cout << "Ingrese un valor aprox. a la raiz (any format): ";

            aprox = this->readNumeric_input(1);

            if (aprox.size() > 0)
            {
                break;
            }

        } while (1);

        do
        {
            std::cout << "Ingrese el error de tolerancia: ";
            err = this->readNumeric_input(1);

            if (err.size() > 0)
                break;
        } while (1);

        metod.NewtonM(aprox[0], err[0], true, true);

        std::cout << "¿Deseas probar otra aproximacion inicial?(S/N) : ";
        if (!this->read_input("s"))
            break;

    } while (1);
}

void Menu::Menu_Secante(int equa)
{
    EquationS metod(4, equa);

    std::vector<float> n, n_1, err;

    do
    {
        std::system(LIMPIAR);
        if (equa == 1)
        {
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "|                    Secante                               |" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "| 1 . f(x) = -SEN(X) + COS(X)                              |" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
        }
        else if (equa == 2)
        {
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "|                    Secante                               |" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "| 2 . f(x) =  e^x*SEN(X)                                   |" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
        }
        else if (equa == 3)
        {
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "|                    Secante                               |" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "| 3 . f(x) =  (x^2 - 5)*COS^3(X)                           |" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
        }
        else if (equa == 4)
        {
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "|                    Secante                               |" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
            std::cout << "| 4 . f(x) =  (x+3)*SEN^2(X)                               |" << std::endl;
            std::cout << "------------------------------------------------------------" << std::endl;
        }
        do
        {
            std::cout << "Ingrese Xn: ";

            n = this->readNumeric_input(1);

            if (n.size() > 0)
            {
                break;
            }

        } while (1);

        do
        {
            std::cout << "Ingrese Xn-1: ";

            n_1 = this->readNumeric_input(1);

            if (n_1.size() > 0)
            {
                break;
            }

        } while (1);

        do
        {
            std::cout << "Ingrese el error de tolerancia: ";
            err = this->readNumeric_input(1);

            if (err.size() > 0)
                break;
        } while (1);

        metod.SecanteM(n[0], n_1[0], err[0], true, true);

        std::cout << "¿Deseas probar otra aproximacion inicial?(S/N) : ";
        if (!this->read_input("s"))
            break;

    } while (1);
}

void Menu::Menu_sistemas_ecuaciones()
{
    do
    {
        std::system(LIMPIAR);
        std::cout << "------------------------------------------------------------" << std::endl;
        std::cout << "|            SISTEMAS DE ECUACIONES LINEALES               |" << std::endl;
        std::cout << "------------------------------------------------------------" << std::endl;
        std::cout << "| 1 . METODOS EXACTOS                                      |" << std::endl;
        std::cout << "|----------------------------------------------------------|" << std::endl;
        std::cout << "| 2 . METODOS ITERATIVOS                                   |" << std::endl;
        std::cout << "------------------------------------------------------------" << std::endl;

        std::vector<float> v;
        do
        {

            std::cout << "Ingrese una opcion del menu(numero): "; // << std::endl;

            v = this->readNumeric_input(1);

            if (v.size() >= 1 && ((((int)v[0]) >= 1) && (((int)v[0]) <= 2)))
            {
                break;
            }
        } while (1);

        this->op_submenu_2 = ((int)v[0]);

        switch (((int)v[0])) //Metodo seleccionado
        {
        case 1:
            this->Menu_metodos_exactos();
            break;

        case 2:
            this->Menu_metodos_iterativos();
            break;

        default:
            break;
        }

        std::cout << "¿Deseas volver al menu anterior?(S/N) : ";
        if (!this->read_input("s"))
            break;

    } while (1);
}

void Menu::Menu_metodos_exactos()
{
    do
    {
        std::system(LIMPIAR);
        std::cout << "------------------------------------------------------------" << std::endl;
        std::cout << "|                SOLUCIONES DE ECUACIONES                  |" << std::endl;
        std::cout << "------------------------------------------------------------" << std::endl;
        std::cout << "| 1 . INVERSION DE MATRICES PARTICIONADO                   |" << std::endl;
        std::cout << "|----------------------------------------------------------|" << std::endl;
        std::cout << "| 2 . GAUSS-JORDAN PARTICIONADO                            |" << std::endl;
        std::cout << "|----------------------------------------------------------|" << std::endl;
        std::cout << "| 3 . METODO DE INTERCAMBIO                                |" << std::endl;
        std::cout << "------------------------------------------------------------" << std::endl;

        std::vector<float> v;
        int op;
        do
        {

            std::cout << "Ingrese una opcion del menu(numero): "; // << std::endl;

            v = this->readNumeric_input(1);

            if (v.size() >= 1 && ((((int)v[0]) >= 1) && (((int)v[0]) <= 3)))
            {
                break;
            }
        } while (1);

        this->op_metodos_Exactos = ((int)v[0]);

        switch (((int)v[0])) //Metodo seleccionado
        {
        case 1:
            this->Menu_inversa_particionado();
            break;
        case 2:
            this->Menu_Gausss_JordanP();
            break;

        default:
            break;
        }

        std::cout << "¿Deseas probar otro metodo?(S/N) : ";
        if (!this->read_input("s"))
            break;

    } while (1);
}

void Menu::Menu_inversa_particionado()
{

    do
    {
        std::vector<std::vector<float>> Mat, vector_indp, t;
        std::vector<float> v;
        float determ;
        int siz;
        std::system(LIMPIAR);
        std::cout << "Ingrese el tamaño de la matriz(un solo numero sin espacios) : ";
        v = this->readNumeric_input(1);
        do
        {
            std::cout << "Ingresa La matriz de coeficientes : " << std::endl;
            try
            {
                Mat = this->read_matrix<float>(((int)v[0]), ((int)v[0]), false);
                break;
            }
            catch (const std::exception &e)
            {
                std::cout << "Invalid Input try again :(" << std::endl;
            }

        } while (1);

        do
        {
            std::cout << "Ingrese el vector de terminos independientes(en forma horizontal): " << std::endl;
            try
            {
                vector_indp = this->read_matrix<float>(1, ((int)v[0]), false);
                break;
            }
            catch (const std::exception &e)
            {
                std::cout << "Invalid Input try again :(" << std::endl;
            }

        } while (1);

        try
        {
            std::cout << "Matriz de Coeficientes " << std::endl;
            print<float>(Mat);
            t = transp(vector_indp);

            std::cout << "Vector de terminos independientes " << std::endl;
            print<float>(t);

            determ = Determinante(Mat);
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
        }

        if (determ == 0)
        {
            std::cout << "La matriz no tiene solucion determinante iguala 0 !" << std::endl;
        }
        else
        {
            std::cout << "Determinante : " << determ << std::endl;
            if (Mat.size() >= 4)
                Inverse_partition(Mat, t, 2, 2);
        }

        std::cout << "¿Deseas probar otra matriz?(S/N) : ";
        if (!this->read_input("s"))
            break;
    } while (1);
}

void Menu::Menu_Gausss_JordanP()
{
    do
    {
        std::vector<std::vector<float>> Mat, vector_indp, t;
        std::vector<float> v;
        float determ;
        int siz;
        std::system(LIMPIAR);
        std::cout << "Ingrese el tamaño de la matriz(un solo numero sin espacios) : ";
        v = this->readNumeric_input(1);
        do
        {
            std::cout << "Ingresa La matriz de coeficientes : " << std::endl;
            try
            {
                Mat = this->read_matrix<float>(((int)v[0]), ((int)v[0]), false);
                break;
            }
            catch (const std::exception &e)
            {
                std::cout << "Invalid Input try again :(" << std::endl;
            }

        } while (1);

        do
        {
            std::cout << "Ingrese el vector de terminos independientes(en forma horizontal): " << std::endl;
            try
            {
                vector_indp = this->read_matrix<float>(1, ((int)v[0]), false);
                break;
            }
            catch (const std::exception &e)
            {
                std::cout << "Invalid Input try again :(" << std::endl;
            }

        } while (1);

        try
        {
            std::cout << "Matriz de Coeficientes " << std::endl;
            print<float>(Mat);
            t = transp(vector_indp);

            std::cout << "Vector de terminos independientes " << std::endl;
            print<float>(t);

            determ = Determinante(Mat);
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
        }

        if (determ == 0)
        {
            std::cout << "La matriz no tiene solucion determinante iguala 0 !" << std::endl;
        }
        else
        {
            std::cout << "Determinante : " << determ << std::endl;
            if (Mat.size() >= 4)
                Gauss_partition(Mat, t, 2, 2);
        }

        std::cout << "¿Deseas probar otra matriz?(S/N) : ";
        if (!this->read_input("s"))
            break;
    } while (1);
}

void Menu::Menu_Metodo_intercambio()
{
}

void Menu::Menu_metodos_iterativos()
{
    do
    {
        std::system(LIMPIAR);
        std::cout << "------------------------------------------------------------" << std::endl;
        std::cout << "|                SISTEMAS DE ECUACIONES                    |" << std::endl;
        std::cout << "------------------------------------------------------------" << std::endl;
        std::cout << "| 1 . METODO DE JACOBI                                     |" << std::endl;
        std::cout << "|----------------------------------------------------------|" << std::endl;
        std::cout << "| 2 . METODO DE GAUSS-SEIDEL                               |" << std::endl;
        std::cout << "|----------------------------------------------------------|" << std::endl;
        std::cout << "| 3 . METODO DE RELAJACION                                 |" << std::endl;
        std::cout << "------------------------------------------------------------" << std::endl;

        std::vector<float> v;
        int op;
        do
        {

            std::cout << "Ingrese una opcion del menu(numero): "; // << std::endl;

            v = this->readNumeric_input(1);

            if (v.size() >= 1 && ((((int)v[0]) >= 1) && (((int)v[0]) <= 3)))
            {
                break;
            }
        } while (1);

        this->op_metodos_Exactos = ((int)v[0]);

        switch (((int)v[0])) //Metodo seleccionado
        {
        case 1:
            this->Menu_jacobi();
            break;

        default:
            break;
        }

        std::cout << "¿Deseas probar otro metodo?(S/N) : ";
        if (!this->read_input("s"))
            break;

    } while (1);
}

void Menu::Menu_jacobi()
{

    do
    {
        std::vector<std::vector<float>> Mat, vector_indp, t;
        std::vector<float> v, error;
        float determ;
        int siz;
        std::system(LIMPIAR);
        std::cout << "Ingrese el tamaño de la matriz(un solo numero sin espacios) : ";
        v = this->readNumeric_input(1);
        do
        {
            std::cout << "Ingresa La matriz de coeficientes : " << std::endl;
            try
            {
                Mat = this->read_matrix<float>(((int)v[0]), ((int)v[0]), false);
                break;
            }
            catch (const std::exception &e)
            {
                std::cout << "Invalid Input try again :(" << std::endl;
            }

        } while (1);

        do
        {
            std::cout << "Ingrese el vector de terminos independientes(en forma horizontal): " << std::endl;
            try
            {
                vector_indp = this->read_matrix<float>(1, ((int)v[0]), false);
                break;
            }
            catch (const std::exception &e)
            {
                std::cout << "Invalid Input try again :(" << std::endl;
            }

        } while (1);

        std::cout << "Ingrese el error de tolerancia" << std::endl;
        error = this->readNumeric_input(1);

        try
        {
            std::cout << "Matriz de Coeficientes " << std::endl;
            print<float>(Mat);
            t = transp(vector_indp);

            std::cout << "Vector de terminos independientes " << std::endl;
            print<float>(t);

            determ = Determinante(Mat);
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
        }

        if (determ == 0)
        {
            std::cout << "La matriz no tiene solucion determinante iguala 0 !" << std::endl;
        }
        else
        {
            std::cout << "Determinante : " << determ << std::endl;
            Jacobi(Mat, t, error[0]);
        }

        std::cout << "¿Deseas probar otra matriz?(S/N) : ";
        if (!this->read_input("s"))
            break;
    } while (1);
}

void Menu::Menu_gaussseidel()
{

    do
    {
        std::vector<std::vector<float>> Mat, vector_indp, t;
        std::vector<float> v, error;
        float determ;
        int siz;
        std::system(LIMPIAR);
        std::cout << "Ingrese el tamaño de la matriz(un solo numero sin espacios) : ";
        v = this->readNumeric_input(1);
        do
        {
            std::cout << "Ingresa La matriz de coeficientes : " << std::endl;
            try
            {
                Mat = this->read_matrix<float>(((int)v[0]), ((int)v[0]), false);
                break;
            }
            catch (const std::exception &e)
            {
                std::cout << "Invalid Input try again :(" << std::endl;
            }

        } while (1);

        do
        {
            std::cout << "Ingrese el vector de terminos independientes(en forma horizontal): " << std::endl;
            try
            {
                vector_indp = this->read_matrix<float>(1, ((int)v[0]), false);
                break;
            }
            catch (const std::exception &e)
            {
                std::cout << "Invalid Input try again :(" << std::endl;
            }

        } while (1);

        std::cout << "Ingrese el error de tolerancia" << std::endl;
        error = this->readNumeric_input(1);

        try
        {
            std::cout << "Matriz de Coeficientes " << std::endl;
            print<float>(Mat);
            t = transp(vector_indp);

            std::cout << "Vector de terminos independientes " << std::endl;
            print<float>(t);

            determ = Determinante(Mat);
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
        }

        if (determ == 0)
        {
            std::cout << "La matriz no tiene solucion determinante iguala 0 !" << std::endl;
        }
        else
        {
            std::cout << "Determinante : " << determ << std::endl;
            GaussSeidel(Mat, t, error[0]);
        }

        std::cout << "¿Deseas probar otra matriz?(S/N) : ";
        if (!this->read_input("s"))
            break;
    } while (1);
}

void Menu::Menu_factorizacionLU()
{
    do
    {
        std::system(LIMPIAR);
        std::cout << "------------------------------------------------------------" << std::endl;
        std::cout << "|            FACTORIZACION LU                              |" << std::endl;
        std::cout << "------------------------------------------------------------" << std::endl;
        std::cout << "| 1 . METODO DE CHOLESKY                                   |" << std::endl;
        std::cout << "|----------------------------------------------------------|" << std::endl;
        std::cout << "| 2 . METODOS DOOLITTLE                                    |" << std::endl;
        std::cout << "|----------------------------------------------------------|" << std::endl;
        std::cout << "| 3 . METODOS CROUT                                        |" << std::endl;
        std::cout << "------------------------------------------------------------" << std::endl;

        std::vector<float> v;
        do
        {

            std::cout << "Ingrese una opcion del menu(numero): "; // << std::endl;

            v = this->readNumeric_input(1);

            if (v.size() >= 1 && ((((int)v[0]) >= 1) && (((int)v[0]) <= 3)))
            {
                break;
            }
        } while (1);

        this->op_submenu_2 = ((int)v[0]);

        switch (((int)v[0])) //Metodo seleccionado
        {
        case 1:
            this->Menu_cholesky();
            break;

        case 2:
            this->Menu_Doolittle();
            break;

        default:
            break;
        }

        std::cout << "¿Deseas volver al menu anterior?(S/N) : ";
        if (!this->read_input("s"))
            break;

    } while (1);
}

void Menu::Menu_cholesky()
{

    do
    {
        std::vector<std::vector<float>> Mat, vector_indp, t;
        std::vector<float> v;
        float determ;
        int siz;
        std::system(LIMPIAR);
        std::cout << "Ingrese el tamaño de la matriz(un solo numero sin espacios) : ";
        v = this->readNumeric_input(1);
        do
        {
            std::cout << "Ingresa La matriz de coeficientes : " << std::endl;
            try
            {
                Mat = this->read_matrix<float>(((int)v[0]), ((int)v[0]), false);
                break;
            }
            catch (const std::exception &e)
            {
                std::cout << "Invalid Input try again :(" << std::endl;
            }

        } while (1);

        do
        {
            std::cout << "Ingrese el vector de terminos independientes(en forma horizontal): " << std::endl;
            try
            {
                vector_indp = this->read_matrix<float>(1, ((int)v[0]), false);
                break;
            }
            catch (const std::exception &e)
            {
                std::cout << "Invalid Input try again :(" << std::endl;
            }

        } while (1);

        try
        {
            std::cout << "Matriz de Coeficientes " << std::endl;
            print<float>(Mat);
            t = transp(vector_indp);

            std::cout << "Vector de terminos independientes " << std::endl;
            print<float>(t);

            //determ = Determinante(Mat);
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
        }

        if (is_simetric(Mat) == false || is_positive(Mat) == false)
        {
            std::cout << "No se puede aplicar el metodo!" << std::endl;
        }
        else
        {
            //std::cout << "Determinante : " << determ << std::endl;
            cholesky<float>(Mat, t);
        }

        std::cout << "¿Deseas probar otra matriz?(S/N) : ";
        if (!this->read_input("s"))
            break;
    } while (1);
}

void Menu::Menu_Doolittle()
{

    do
    {
        std::vector<std::vector<float>> Mat, vector_indp, t;
        std::vector<float> v;
        float determ;
        int siz;
        std::system(LIMPIAR);
        std::cout << "Ingrese el tamaño de la matriz(un solo numero sin espacios) : ";
        v = this->readNumeric_input(1);
        do
        {
            std::cout << "Ingresa La matriz de coeficientes : " << std::endl;
            try
            {
                Mat = this->read_matrix<float>(((int)v[0]), ((int)v[0]), false);
                break;
            }
            catch (const std::exception &e)
            {
                std::cout << "Invalid Input try again :(" << std::endl;
            }

        } while (1);

        do
        {
            std::cout << "Ingrese el vector de terminos independientes(en forma horizontal): " << std::endl;
            try
            {
                vector_indp = this->read_matrix<float>(1, ((int)v[0]), false);
                break;
            }
            catch (const std::exception &e)
            {
                std::cout << "Invalid Input try again :(" << std::endl;
            }

        } while (1);

        try
        {
            std::cout << "Matriz de Coeficientes " << std::endl;
            print<float>(Mat);
            t = transp(vector_indp);

            std::cout << "Vector de terminos independientes " << std::endl;
            print<float>(t);

            //determ = Determinante(Mat);
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << '\n';
        }

        //std::cout << "Determinante : " << determ << std::endl;
        Doolittle<float>(Mat, t);

        std::cout << "¿Deseas probar otra matriz?(S/N) : ";
        if (!this->read_input("s"))
            break;
    } while (1);
}

Menu::~Menu()
{
}

int main(int argc, char const *argv[])
{

    Menu m;

    while (1)
        if (m.Main_menu() == 4)
            break;

    return 0;
}
