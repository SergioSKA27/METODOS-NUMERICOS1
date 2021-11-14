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

    std::vector<float> readNumeric_input(int argsN); //lee una linea y la retorna como un string
    bool read_input(std::string find);               // lee un input y busca una subcadena(find) retorna true si existe una ocurrencia

public:
    Menu(/* args */);

    int Main_menu();

    void Menu_Soluciones_de_ecuaciones();
    int Menu_ecuaciones();
    void Menu_Biseccion(int equa);
    void Menu_FalsaPos(int equa);
    void Menu_Newton(int equa);
    void Menu_Secante(int equa);

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
    //*si el usuario introduce mas de los N argumentos dichos seran descartados
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
