#include <bits/stdc++.h>
#ifdef __linux__
#define LIMPIAR "clear"
#endif // __linux__

#ifdef __MINGW32__
#define LIMPIAR "CLS"
#endif // __MINGW32__

class EquationS
{

private:
    int equation;
    int method;

    inline float evaluar(float value, int func, bool derivate); // Evalua una funcion de acuerdo a un valor
    void bisection_method(std::pair<float, float> &interval, float error, float &root, bool print, int k);

public:
    EquationS(int method, int eq);

    float BiseccionM(std::pair<float, float> interval, float error, bool print);

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

void EquationS::bisection_method(std::pair<float, float> &interval, float error, float &root, bool print, int k)
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

    this->bisection_method(ninter, error, root, print, k + 1);
    return;
}

float EquationS::BiseccionM(std::pair<float, float> interval, float error, bool print)
{
    float fb = this->evaluar(interval.second, this->equation, false); // f(b)
    float fa = this->evaluar(interval.first, this->equation, false);  // f(a)
    float dfb = this->evaluar(interval.second, this->equation, true); // f'(b)
    float dfa = this->evaluar(interval.first, this->equation, true);  // f'(a)
    float root = 0.0;

    std::cout << "f(" << interval.first << ") = " << fa << ", f(" << interval.second << ") = " << fb << std::endl;
    std::cout << "f'(" << interval.first << ") = " << dfa << ", f'(" << interval.second << ") = " << dfb << std::endl;

    if (!(fa * fb < 0 && dfa * dfb > 0))
        return 0;
    //throw std::invalid_argument(":("); //No hay raiz aislada

    std::cout << "|  a   |  b  |   f(a)   |   f(b)   | x = (a+b)/2 |   |f(x)| < error   |   f(x)*f(a) < 0   | " << std::endl;

    this->bisection_method(interval, error, root, print, 0);

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

public:
    Menu(/* args */);

    void Main_menu();
    void Menu_Soluciones_de_ecuaciones();
    void Menu_ecuaciones();

    ~Menu();
};

Menu::Menu()
{
    this->op_main_menu = 0;
    this->op_sub_menu1 = 0;
    this->function_submenu1 = 0;
}

void Menu::Main_menu()
{
    std::system(LIMPIAR);
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << "|                   METODOS NUMERICOS                      |" << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << "| 1 . SOLUCIONES DE ECUACIONES                             |" << std::endl;
    std::cout << "|----------------------------------------------------------|" << std::endl;
    std::cout << "| 2 . SISTEMAS DE ECUACIONES LINEALES                      |" << std::endl;
    std::cout << "|----------------------------------------------------------|" << std::endl;
    std::cout << "| 3 . FACTORIZACION LU                                     |" << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;

    std::cout << "Ingrese una opcion del menu(numero): "; // << std::endl;

    while (std::cin >> this->op_main_menu && (this->op_main_menu < 1 || this->op_main_menu > 3))
    {
        std::cout << "Ingrese una Opcion del menu: ";
    }

    std::system(LIMPIAR);

    switch (this->op_main_menu)
    {
    case 1:

        this->Menu_Soluciones_de_ecuaciones();

        break;

    default:
        break;
    }
}

void Menu::Menu_Soluciones_de_ecuaciones()
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

    std::cout << "Ingrese una opcion del menu(numero): "; // << std::endl;

    while (std::cin >> this->op_sub_menu1 && (this->op_sub_menu1 < 1 || this->op_sub_menu1 > 4))
    {
        std::cout << "Ingrese una Opcion del menu: ";
    }

    if (this->op_sub_menu1 == 1) //metodo de biseccion
    {
        std::string ex;
        do
        {
            int eq;
            std::system(LIMPIAR);

            std::string rep;

            this->Menu_ecuaciones();
            while (std::cin >> eq && (eq < 1 || eq > 4))
            {
                std::cout << "Ingrese una Opcion del menu: ";
            }

            do
            {
                std::string s1 = "", s2 = "";
                std::pair<float, float> interval;
                char s[100];
                int a, b;
                int k = 0;
                float error;
                std::string im, p;
                bool print;
                std::system(LIMPIAR);

                if (eq == 1)
                {
                    std::cout << "------------------------------------------------------------" << std::endl;
                    std::cout << "| 1 . f(x) = -SEN(X) + COS(X)                              |" << std::endl;
                    std::cout << "------------------------------------------------------------" << std::endl;
                }
                else if (eq == 2)
                {
                    std::cout << "------------------------------------------------------------" << std::endl;
                    std::cout << "| 2 . f(x) =  e^x*SEN(X)                                   |" << std::endl;
                    std::cout << "------------------------------------------------------------" << std::endl;
                }
                else if (eq == 3)
                {
                    std::cout << "------------------------------------------------------------" << std::endl;
                    std::cout << "| 3 . f(x) =  (x^2 - 5)*COS^3(X)                           |" << std::endl;
                    std::cout << "------------------------------------------------------------" << std::endl;
                }
                else if (eq == 4)
                {
                    std::cout << "------------------------------------------------------------" << std::endl;
                    std::cout << "| 4 . f(x) =  (x+3)*SEN^2(X)                               |" << std::endl;
                    std::cout << "------------------------------------------------------------" << std::endl;
                }

                std::cout << "INGRESE EL INTERVALO(any format):  ";
                std::cin.ignore();
                std::cin.getline(s, 100, '\n');
                std::cout << "INGRESE EL ERROR DE TOLERANCIA : ";
                std::cin >> error;
                std::cout << "Desea Imprimir el metodo en pantalla(S/N): ";
                std::cin.ignore();
                std::cin >> im;

                std::transform(im.begin(), im.end(), im.begin(), tolower);
                std::cout << s << std::endl;
                if (*im.begin() == 's')
                    print = true;
                else
                {
                    print = false;
                }

                for (size_t i = 0; s[i]; i++)
                {
                    if ((s[i] - '0' >= 0 && s[i] - '0' <= 9) || s[i] == '.')
                        s1 += s[i];

                    if (s[i] == ',' || s[i] == ' ')
                    {
                        k = i;
                        //std::cout << i << endl;
                        break;
                    }
                }

                for (size_t i = k; s[i]; i++)
                {
                    if ((s[i] - '0' >= 0 && s[i] - '0' <= 9) || s[i] == '.')
                        s2 += s[i];
                }

                std::cout << s1 << " -<- " << s2 << std::endl;
                interval.first = std::stof(s1);
                interval.second = std::stof(s2);
                EquationS bis(1, eq);

                bis.BiseccionM(interval, error, print);

                std::cout << "Deseas evaluar otro intervalo(S/N): ";
                std::cin.ignore();
                std::getline(std::cin, rep);
                std::transform(rep.begin(), rep.end(), rep.begin(), tolower);

            } while (*rep.begin() == 's' || *rep.begin() == 'y');

            std::cout << "Deseas regresar al menu anterior(S/N): ";
            std::cin.ignore();
            std::getline(std::cin, ex);
            std::transform(ex.begin(), ex.end(), ex.begin(), tolower);

        } while (*ex.begin() == 's' || *ex.begin() == 'y');
    }
}

void Menu::Menu_ecuaciones()
{
    std::system(LIMPIAR);
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << "|                ECUACIONES                                |" << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << "| 1 . f(x) = -SEN(X) + COS(X)                              |" << std::endl;
    std::cout << "|----------------------------------------------------------|" << std::endl;
    std::cout << "| 2 . f(x) =  e^x*SEN(X)                                   |" << std::endl;
    std::cout << "|----------------------------------------------------------|" << std::endl;
    std::cout << "| 3 . f(x) =  (x^2 - 5)*COS^3(X)                           |" << std::endl;
    std::cout << "|----------------------------------------------------------|" << std::endl;
    std::cout << "| 4 . f(x) =  (x+3)*SEN^2(X)                               |" << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;
}

Menu::~Menu()
{
}

int main(int argc, char const *argv[])
{

    Menu m;
    while (1)
    {
        m.Main_menu();
    }

    return 0;
}
