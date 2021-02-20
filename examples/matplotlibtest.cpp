

#include <matplotlibcpp.h>
#include <vector>

int main(int argc, char *argv[]) {
    std::cout << "HELLP\n";
    std::vector<double> x;
    std::vector<double> f;

    for (int i = 0; i < 100; i++) {
        x.push_back(3.2*i/100.0);
        f.push_back(std::sin(x[i]));
    }

    matplotlibcpp::plot(x, f);

    matplotlibcpp::show();
}