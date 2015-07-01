#include <iostream>

int main()
{
    #ifdef __ICC
    std::cout << "using icc" << std::endl;
    std::cout << __VERSION__ << std::endl;
    #else
    std::cout << "using gcc" << std::endl;
    #endif

    return 0;
}
