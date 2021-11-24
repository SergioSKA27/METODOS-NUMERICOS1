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

/*
    Creamos una clase Numero para poder efectuar un mayor numero de opraciones aritmeticas 
    ya sea con numeros reales, enteros , racionales e irracionales.
    Ademas dicha clase nos permite expresar el resultado de una opracion de la manera exacta 
    haciendo uso de fracciones para expresar los resultados.
*/

class Number
{
private:
    float real_part;
    int numerator, denominator; //Representacion como racional
    bool is_rational;
    bool is_constant;

    /*
        podemos usar alguna de las constantes mas conocidas 
        declarando el objeto con su respectivo nombre:
            
            - π  = PI
            - e = E o e
            - √2  = SQRT2
            - Φ = PHI
    
        *Nota: las claves unicamente deben contener dichos caracteres ya sea mayusculas o minusculas
    */

public:
    Number();
    Number(float real);       //Constructor reales
    Number(float a, float b); //Constructor racional
    Number(std::string key);  //constructor Irracional
    //fracciones con irracionales
    Number(std::string numerator, float denominator);       //  Irracional/real
    Number(std::string numerator, std::string denominator); //  Irracional/Irracional
    Number(float numerator, std::string denominator);       //  real/Irracional

    friend std::ostream &operator<<(std::ostream &o, const Number n)
    {
        if (n.is_rational)
        {
            if (n.numerator > 0 && n.denominator > 0)
                o << n.numerator << "/" << n.denominator;
            else if (n.numerator < 0 && n.denominator < 0)
                o << -1 * n.numerator << "/" << -1 * n.denominator;
            else
                o << '-' << abs(n.numerator) << "/" << abs(n.denominator);
        }
        else
            o << n.real_part;

        return o;
    }

    friend std::ostream &operator<<(std::ostream &o, const Number *n)
    {
        if (n->is_rational)
        {
            if (n->numerator > 0 && n->denominator > 0)
                o << n->numerator << "/" << n->denominator;
            else if (n->numerator < 0 && n->denominator < 0)
                o << -1 * n->numerator << "/" << -1 * n->denominator;
            else
                o << '-' << abs(n->numerator) << "/" << abs(n->denominator);
        }
        else
            o << n->real_part;

        return o;
    }

    static void *operator new(size_t size);
    static void *operator new[](size_t size);
    void operator delete[](void *p);
    void operator delete(void *p);

    Number &simplify();

    Number &operator=(const Number &x);
    Number &operator=(const float x);
    Number &operator=(const std::pair<float, float> &x);

    Number operator+(const Number &x);
    Number operator+(const float &x);

    Number &operator+=(const Number &x);
    Number &operator+=(const float &x);

    Number operator-(const Number &x);
    Number operator-(const float &x);

    Number operator*(const Number &x);
    Number operator*(const float &x);

    Number operator/(const Number &x);
    Number operator/(const float &x);

    Number make_fraction(float &a, float &b);
    ~Number();
};

Number::Number()
{
    this->is_constant = false;
    this->is_rational = false;
    this->real_part = 0;
    this->numerator = 0;
    this->denominator = 1;
}

Number::Number(float real)
{
    this->is_constant = false;
    this->is_rational = false;
    this->real_part = real;
    this->numerator = real;
    this->denominator = 1;
}

Number::Number(float a, float b)
{
    if (b == 0)
        throw std::invalid_argument("Division por 0 !");
    this->is_constant = false;
    if (b != 1)
        this->is_rational = true;
    this->real_part = a / b;
    this->numerator = a;
    this->denominator = b;
}

Number::Number(std::string key)
{
    const float Pi = M_PI, E = exp(1);
    const float sqrt2 = M_SQRT2;
    size_t pos, k;
    std::string claves[] = {"pi", "e", "sqrt2", "phi"};

    this->is_constant = true;
    this->is_rational = false;
    std::transform(key.begin(), key.end(), key.begin(), tolower);

    for (size_t i = 0; i < 4; i++)
    {
        pos = key.find(claves[i]);
        if (pos != std::string::npos)
        {
            k = i;
            break;
        }
    }

    if (pos == std::string::npos)
        throw std::invalid_argument((key + " No existe :("));
    else
    {

        if (k == 0)
        {
            this->real_part = Pi;
            this->numerator = Pi;
        }

        if (k == 1)
        {
            this->real_part = E;
            this->numerator = E;
        }

        if (k == 3)
        {
            this->real_part = sqrt2;
            this->numerator = sqrt2;
        }
        if (k == 4)
        {
            this->real_part = 1.61803398875;
            this->numerator = 1.61803398875;
        }
    }

    this->denominator = 1;
}

Number &Number::simplify()
{
    bool is_simplifying = false;
    if (this->numerator == this->denominator)
    {
        this->numerator = 1;
        this->denominator = 1;
        this->is_rational = false;
        return *this;
    }
    else
    {
        for (int i = 2; i < this->numerator + 1; i++)
        {
            if (this->numerator % i == 0)
            {
                if (this->denominator % i == 0)
                {
                    this->numerator = this->numerator / i;
                    this->denominator = this->denominator / i;
                    is_simplifying = true;
                    break;
                }
            }
        }
        if (!is_simplifying)
            return *this;
        else
            this->simplify();
    }
    return *this;
}

Number &Number::operator=(const Number &x)
{
    this->denominator = x.denominator;
    this->numerator = x.numerator;
    this->is_constant = x.is_constant;
    this->is_rational = x.is_rational;
    this->real_part = x.real_part;

    return *this;
}

Number &Number::operator=(const float x)
{
    this->denominator = 1;
    this->numerator = x;
    this->is_constant = false;
    this->is_rational = false;
    this->real_part = x;

    return *this;
}

Number Number::operator+(const Number &x)
{

    if (this->is_rational || x.is_rational)
    {
        if (x.denominator == this->denominator)
        {
            Number r(this->numerator + x.numerator, this->denominator);
            return r.simplify();
        }
        else
        {
            float a = this->numerator * x.denominator;
            float b = this->denominator * x.numerator;
            Number r(a + b, this->denominator * x.denominator);
            return r.simplify();
        }
    }

    float y = this->real_part + x.real_part;
    Number r(y);

    return r;
}

Number Number::operator+(const float &x)
{

    if (this->is_rational || (((int)x) == x))
    {
        if (this->denominator == 1)
        {
            Number r(this->numerator + x, this->denominator);
            return r.simplify();
        }
        else
        {
            float a = this->numerator;
            float b = this->denominator * x;
            Number r(a + b, this->denominator * 1);
            return r.simplify();
        }
    }

    float y = this->real_part + x;
    Number r(y);

    return r;
}

Number &Number::operator+=(const Number &x)
{

    if (this->is_rational || x.is_rational)
    {
        if (x.denominator == this->denominator)
        {
            this->numerator += x.numerator;

            return *this;
        }
        else
        {
            float a = this->numerator * x.denominator;
            float b = this->denominator * x.numerator;
            this->numerator = a + b;
            this->denominator = this->denominator * x.denominator;
            return *this;
        }
    }

    float y = this->real_part + x.real_part;
    Number r(y);

    (*this) = r;

    return *this;
}

Number &Number::operator+=(const float &x)
{

    if (this->is_rational || (((int)x) == x))
    {
        if (this->denominator == 1)
        {
            Number r(this->numerator + x, this->denominator);
            r.simplify();
            (*this) = r;
            return *this;
        }
        else
        {
            float a = this->numerator;
            float b = this->denominator * x;
            Number r(a + b, this->denominator * 1);
            r.simplify();
            (*this) = r;
            return *this;
        }
    }

    float y = this->real_part + x;
    Number r(y);
    (*this) = r;
    return *this;
}

Number Number::operator-(const Number &x)
{

    if (this->is_rational || x.is_rational)
    {
        if (x.denominator == this->denominator)
        {
            Number r(this->numerator - x.numerator, this->denominator);
            return r.simplify();
        }
        else
        {
            float a = this->numerator * x.denominator;
            float b = this->denominator * x.numerator;
            Number r(a - b, this->denominator * x.denominator);
            return r.simplify();
        }
    }

    float y = this->real_part - x.real_part;
    Number r(y);

    return r;
}

Number Number::operator-(const float &x)
{

    if (this->is_rational || (((int)x) == x))
    {
        if (this->denominator == 1)
        {
            Number r(this->numerator - x, this->denominator);
            return r.real_part;
        }
        else
        {
            float a = this->numerator;
            float b = this->denominator * x;
            Number r(a - b, this->denominator * 1);
            return r.real_part;
        }
    }

    float y = this->real_part - x;
    Number r(y);
    return r;
}

Number Number::operator*(const Number &x)
{

    if (this->is_rational || x.is_rational)
    {
        Number r(this->numerator * x.numerator, this->denominator * x.denominator);
        return r.simplify();
    }

    float y = this->real_part * x.real_part;
    Number r(y);

    return r;
}

Number Number::operator*(const float &x)
{

    if (this->is_rational || (((int)x) == x))
    {

        Number r(this->numerator * x, this->denominator);
        return r.real_part;
    }

    float y = this->real_part * x;
    Number r(y);

    return r;
}

Number Number::operator/(const Number &x)
{
    if (x.real_part == 0)
    {
        throw std::invalid_argument("Division por 0!");
    }

    if (this->is_rational || x.is_rational)
    {

        float a = this->denominator * x.numerator;
        float b = this->numerator * x.denominator;
        Number r(a, b);
        return r.simplify();
    }

    float y = this->real_part / x.real_part;
    Number r(y);

    return r;
}

Number Number::operator/(const float &x)
{
    if (x == 0)
    {
        throw std::invalid_argument("Division por 0!");
    }

    if (this->is_rational || (((int)x) == x))
    {
        float a = this->denominator * x;
        float b = this->numerator;
        Number r(a, b);
        return r.real_part;
    }

    float y = this->real_part / x;
    Number r(y);

    return r;
}

Number &Number::operator=(const std::pair<float, float> &x)
{
    Number r(x.first, x.second);

    (*this) = r;
    return *this;
}

Number Number::make_fraction(float &a, float &b)
{
    Number x(a, b);
    return x;
}

static void *Number::operator new(size_t size)
{
    std::cout << "New" << std::endl;
    void *p = std::malloc(size);

    if (!p)
    {
        std::cerr << "Error al reservar memoria!";
        throw std::bad_alloc();
    }

    std::cout << "New works" << std::endl;

    return p;
}

void Number::operator delete(void *p)
{
    //std::cout << "Delete" << std::endl;
    free(p);
}

static void *Number::operator new[](size_t size)
{
    void *p;

    p = std::malloc(size);

    if (!p)
    {
        throw std::bad_alloc();
    }
    //std::cout << "New works" << std::endl;

    return p;
}

void Number::operator delete[](void *p)
{
    free(p);
}

Number::~Number()
{
}

/*
La clase Matriz funciona como un template es decir que podemos instanciar el objeto 
Matriz con cualquier tipo las operaciones con matrices solo estan disponibles con tipos
Numericos o objetos con los operadores aritmeticos sobrecargados
*/

template <class type>
class Matriz
{
private:
    type **Mat;
    int filas;
    int columnas;

    inline Matriz<type> ExtractMat(type **Mat, int sz, int F, int C); //funcion para el Determinate
    type Det(type **Mat, int sz);                                     //Calcula el determinante por cofactores

    type brackethelp(type *fila, int col); //funcion de ayuda para el operador corche
    int brcketaux, bracketaux2;            //variables auxiliares para indices en corchetes

public:
    Matriz(const type init, int filas, int columnas); //Para cualquier Matriz inicializada
    Matriz(int filas, int columnas);                  //Para cualquier Matriz
    Matriz(int size);                                 //Para Matrices cuadradas
    Matriz();                                         //Constructor vacio

    Matriz &resize(int filas, int columnas); //Redimensionar Matriz sin perder los datos

    void print(); //Imprimir la matriz

    type Determinante(); //Devuelve el determinante de la matriz

    Matriz<type> transp(); //Devuelve la transpuesta de la matriz
    Matriz<type> adj();    //devuelve la adjunta de la matriz

    Matriz<type> inversa(); //Retorna la inversa de una matriz (mediante la ajunta de la matriz)

    Matriz<type> operator+(const Matriz<type> &Mat2);  //Suma de matrices
    Matriz<type> operator-(const Matriz<type> &Mat2);  //Resta de matrices
    Matriz<type> operator*(const Matriz<type> &Mat2);  //Multiplicacion de matrices
    Matriz<type> &operator=(const Matriz<type> &Mat2); //Operador de asignacion

    type *operator[](const int index);        //Operador corchete para filas
    type operator[](short int index2);        //Operador corchete para columnas
    Matriz<type> &operator=(const type Data); //Para asignar valor a las casillas de la matriz

    ~Matriz();
};

template <class type>
Matriz<type>::Matriz(int filas, int columnas)
{
    if (this->filas != 1 && this->columnas != 1)
    {
        this->Mat = new type *[filas];

        for (int i = 0; i < filas; i++)
        {
            this->Mat[i] = new type[columnas];
        }
    }
    this->filas = filas;
    this->columnas = columnas;
}

template <class type>
Matriz<type>::Matriz(int size)
{
    this->Mat = new type *[size];

    for (int i = 0; i < size; i++)
    {
        this->Mat[i] = new type[size];
    }
    this->filas = size;
    this->columnas = size;
}

template <class type>
Matriz<type>::Matriz(const type init, int filas, int columnas)
{
    this->Mat = new type *[filas];

    for (int i = 0; i < filas; i++)
    {
        this->Mat[i] = new type[columnas];
    }

    for (int i = 0; i < filas; i++)
    {
        for (int j = 0; j < columnas; j++)
        {
            this->Mat[i][j] = init;
        }
    }

    this->filas = filas;
    this->columnas = columnas;
}

template <class type>
Matriz<type>::Matriz()
{
    this->Mat = NULL;
    this->filas = 0;
    this->columnas = 0;
}

template <class type>
Matriz<type> &Matriz<type>::resize(int filas, int columnas)
{
    type **newMat;

    newMat = new type *[filas];

    for (int i = 0; i < filas; i++)
    {
        newMat[i] = new type[columnas];
    }

    for (int i = 0; i < this->filas; i++)
    {
        for (int j = 0; j < this->columnas; j++)
        {
            if (i < this->filas && j < this->columnas)
                newMat[i][j] = this->Mat[i][j];
        }
    }

    this->Mat = newMat;
    this->filas = filas;
    this->columnas = columnas;

    return *this;
}

template <class type>
inline Matriz<type> Matriz<type>::ExtractMat(type **Mat, int sz, int F, int C)
{
    Matriz<type> result(sz - 1, sz - 1);
    int k = 0, l = 0;

    for (int i = 0; i < sz; i++)
    {
        l = 0;
        for (int j = 0; j < sz; j++)
        {
            if (i != F && j != C)
            { //Si no estamos en la fila y la columna que se van a eliminar asignamos
                //el valor en esa posicion al valor k,l de la matriz resultado
                result.Mat[k][l] = Mat[i][j];
                l++; //iteramos las columnas de la matriz resultado
            }
        }
        if (i != F) //si no estamos en la fila que se va a eliminar iteramos k
            k++;    //iteramos las filas de la matriz resultado
    }

    return result;
}

template <class type>
type Matriz<type>::Det(type **Mat, int sz)
{

    type detval;
    type dt;
    std::vector<type> v;

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
        Matriz<type> res(sz - 1, sz - 1);

        for (int i = 0; i < sz; i++)
        {
            res = this->ExtractMat(Mat, sz, 0, i);
            // res.print();

            dt = this->Det(res.Mat, sz - 1);

            //std::cout << "det -> " << dt << std::endl;

            if (i % 2 == 0)
            {
                v.push_back(Mat[0][i] * dt);
                //std::cout << Mat[0][i] * dt << std::endl;
                //std::cout << "+" << Mat[0][i] << "*" << dt << std::endl;
            }
            else
            {
                type men;
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

template <class type>
type Matriz<type>::Determinante()
{
    if (this->filas != this->columnas)
        throw std::invalid_argument("La matriz no es cuadrada");

    return this->Det(this->Mat, this->filas);
}

template <class type>
Matriz<type> Matriz<type>::transp()
{
    Matriz<type> T(this->columnas, this->filas);

    for (int i = 0; i < this->columnas; i++)
    {
        for (int j = 0; j < this->filas; j++)
        {
            T.Mat[i][j] = this->Mat[j][i];
        }
    }

    return T;
}

template <class type>
Matriz<type> Matriz<type>::adj()
{
    Matriz<type> AD(this->filas, this->columnas);
    Matriz<type> aux(this->filas - 1, this->columnas - 1);

    for (int i = 0; i < this->columnas; i++)
    {

        for (int j = 0; j < this->filas; j++)
        {
            aux = this->ExtractMat(this->Mat, this->filas, i, j);
            if (((i + 1) + (j + 1)) % 2 == 0)
            {
                AD.Mat[i][j] = aux.Determinante();
            }
            else
            {
                AD.Mat[i][j] = aux.Determinante() * -1;
            }
        }
    }

    return AD.transp();
}

template <class type>
Matriz<type> Matriz<type>::inversa()
{ /*Calculamos la inversa sabiendo que la inversa de una matriz A es:
            trans(adj(A))
    inv(A)= -------------
                |A|
    Donde |A| representa el determiante de la matriz.
    (la  Inversa  de  la  matriz  solo es  aplicable  apartir  del  conjuto 
    de los numeros reales(Naturales, Enteros, Racionales e Irracionales no))
*/
    if (this->Determinante() == 0)
        throw std::invalid_argument("La matriz no tiene inversa(Determinante igual a 0\n");

    Matriz<type> Inv, aux;

    type det = this->Determinante();

    aux = this->adj();

    Inv = aux.transp();

    for (int i = 0; i < Inv.filas; i++)
    {
        for (int j = 0; j < Inv.columnas; j++)
        {
            Inv.Mat[i][j] = Inv.Mat[i][j] / det;
        }
    }

    return Inv;
}

template <class type>
void Matriz<type>::print()
{
    if (Mat == NULL)
        throw std::invalid_argument("La Matriz es de tamano 0\n");

    for (int i = 0; i < this->filas; i++)
    {
        for (int j = 0; j < this->columnas; j++)
        {
            std::cout << this->Mat[i][j] << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

template <class type>
Matriz<type> Matriz<type>::operator+(const Matriz<type> &Mat2)
{
    if (this->filas != Mat2.filas || this->columnas != Mat2.columnas)
        throw std::invalid_argument("Las Matrices no tienen el mismo tamano\n");
    Matriz<type> Result(this->filas, this->columnas);

    for (int i = 0; i < this->filas; i++)
    {
        for (int j = 0; j < this->columnas; j++)
        {
            Result.Mat[i][j] = this->Mat[i][j] + Mat2.Mat[i][j];
        }
    }

    return Result;
}

template <class type>
Matriz<type> Matriz<type>::operator-(const Matriz<type> &Mat2)
{
    if (this->filas != Mat2.filas || this->columnas != Mat2.columnas)
        throw std::invalid_argument("Las Matrices no tienen el mismo tamano\n");
    Matriz<type> Result(this->filas, this->columnas);

    for (int i = 0; i < this->filas; i++)
    {
        for (int j = 0; j < this->columnas; j++)
        {
            Result.Mat[i][j] = this->Mat[i][j] - Mat2.Mat[i][j];
        }
    }

    return Result;
}

template <class type>
Matriz<type> Matriz<type>::operator*(const Matriz<type> &Mat2)
{
    if (this->filas != Mat2.columnas)
        throw std::invalid_argument("Las Matrices no se pueden multiplicar\n");
    Matriz<type> Result(this->filas, Mat2.columnas);

    type sum; //si se planea usar objetos estos tienen que tener un 0 y sobrecargar el operador = para poder asignarlo

    sum = 0;

    for (int i = 0; i < this->filas; i++)
    {
        sum = 0;
        for (int j = 0; j < Mat2.columnas; j++)
        {
            for (int k = 0; k < Mat2.columnas; k++)
            {
                sum = sum + (this->Mat[i][k] * Mat2.Mat[k][j]);
                Result.Mat[k][j] = sum;
            }
        }
    }

    return Result;
}

template <class type>
Matriz<type> &Matriz<type>::operator=(const Matriz<type> &Mat2)
{
    if (this->filas != Mat2.filas || this->columnas != Mat2.columnas)
        this->resize(Mat2.filas, Mat2.columnas);

    for (int i = 0; i < this->filas; i++)
    {
        for (int j = 0; j < this->columnas; j++)
        {
            this->Mat[i][j] = Mat2.Mat[i][j];
        }
    }

    this->filas = Mat2.filas;
    this->columnas = Mat2.columnas;

    return *this;
}

template <class type>
type Matriz<type>::brackethelp(type *fila, int col)
{
    return fila[col];
}

template <class type>
type *Matriz<type>::operator[](const int index)
{
    if (index >= this->filas)
        throw std::out_of_range("Fuera del rango de la matriz\n");
    this->brcketaux = index;
    return this->Mat[index];
}

template <class type>
type Matriz<type>::operator[](short int index2)
{
    if (index2 >= this->columnas)
        throw std::out_of_range("Fuera del rango de la matriz\n");

    this->bracketaux2 = index2;
    return (this->brackethelp((*this)[this->brcketaux], index2));
}

template <class type>
Matriz<type> &Matriz<type>::operator=(const type Data)
{
    this->Mat[brcketaux][bracketaux2] = Data;
    return *this;
}

template <class type>
Matriz<type>::~Matriz()
{
    /*for (int i = 0; i < this->filas; i++)
    {
        delete this->Mat[i];
    }
    delete this->Mat;*/
}

template <class T>
class Linear_equation
{
private:
    Matriz<T> Coef_Matriz;
    Matriz<T> vector_indp;

public:
    Linear_equation(Matriz<T> MatCoef, Matriz<T> vector_indp);
    ~Linear_equation();
};
template <class T>
Linear_equation<T>::Linear_equation(Matriz<T> MatCoef, Matriz<T> vector_indp)
{
    this->Coef_Matriz = MatCoef;
    this->vector_indp = vector_indp.transp();
}
template <class T>
Linear_equation<T>::~Linear_equation()
{
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

    std::vector<float> readNumeric_input(int argsN); //lee una linea y la retorna como un string
    bool read_input(std::string find);               // lee un input y busca una subcadena(find) retorna true si existe una ocurrencia

public:
    Menu();
    template <class F>
    Matriz<F> read_matrix(size_t filas, size_t columnas, bool n); //lee una matriz de n * m y retorna la matriz

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

    void Menu_metodos_iterativos();

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
        if ((line[i] - '0' >= 0 && line[i] - '0' <= 9) || (line[i] == '.'))
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
Matriz<F> Menu::read_matrix(size_t filas, size_t columnas, bool n)
{
    std::vector<std::string> v;
    std::vector<std::vector<float>> m;
    Matriz<F> Mat(filas, columnas);

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
                        Mat[i][k] = std::stof(aux);
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

                        Number frac;
                        if (n)
                        {
                            frac = std::make_pair(nume, den);
                            //Mat[i][k] = frac;
                        }
                        else
                            Mat[i][k] = nume / den;

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

    std::cout << "matrix :" << std::endl;

    Mat.print();

    return Mat;
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

        do
        {
            this->op_submenu_2 = ((int)v[0]);

            switch (((int)v[0])) //Metodo seleccionado
            {
            case 1:
                this->Menu_metodos_exactos();
                break;

            case 2:
                //this->Menu_metodos_iterativos();
                break;

            default:
                break;
            }

            std::cout << "¿Deseas volver al menu anterior?(S/N) : ";
            if (!this->read_input("s"))
                break;
        } while (1);

        std::cout << "¿Deseas probar otro metodo?(S/N) : ";
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

        do
        {
            this->op_metodos_Exactos = ((int)v[0]);

            switch (((int)v[0])) //Metodo seleccionado
            {
            case 0:
                break;

            default:
                break;
            }

            std::cout << "¿Deseas probar otro metodo?(S/N) : ";
            if (!this->read_input("s"))
                break;
        } while (1);

        std::cout << "¿Deseas volver al menu anterior?(S/N) : ";
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

    //while (1)
    //    if (m.Main_menu() == 4)
    //        break;
    Matriz<float> A(2);
    A = m.read_matrix<float>(2, 2, false);

    //A.adj().print();

    std::cout << A.Determinante() << std::endl;

    return 0;
}
