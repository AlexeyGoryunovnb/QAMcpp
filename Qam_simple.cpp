#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <bitset>
#include <fstream>
#include <cstdlib>


double db_to_linear(double db) {
    return std::pow(10.0, db / 10.0);
}

std::vector<int> int_to_bits(int value, int num_bits) {
    std::vector<int> bits(num_bits);
    for (int i = 0; i < num_bits; ++i) {
        bits[i] = (value >> i) & 1;
    }
    return bits;
}

class RandomGenerator {
private:
    std::mt19937 gen;
    std::uniform_int_distribution<> int_dist;
    std::normal_distribution<> normal_dist;

public:
    RandomGenerator() : gen(std::random_device{}()),
        int_dist(0, std::numeric_limits<int>::max()),
        normal_dist(0.0, 1.0) {}

    int rand_int(int min, int max) {
        return int_dist(gen) % (max - min + 1) + min;
    }

    double rand_normal() {
        return normal_dist(gen);
    }
};

class QAMModulator_simple {
private:
    int p;
    int q;
    int A;

public:
    QAMModulator_simple(int p_val)
        : p(p_val) {
        q = std::pow(2 * (std::pow(2, p)), 2);
        A = 2 * std::pow(2, p) - 1;
    }

    int get_q() {
        return q;
    }

    void modulate(int n, double& x, double& y) {
        n = n % q;
        int x_idx = n % 2;
        int y_idx = (n / 2) % 2;

        if (p > 0) {
            x_idx = std::abs(3 * ((n / 4) % 2) - x_idx);
            y_idx = std::abs(3 * ((n / 8) % 2) - y_idx);
            if (p > 1) {
                x_idx = std::abs(7 * ((n / 16) % 2) - x_idx);
                y_idx = std::abs(7 * ((n / 32) % 2) - y_idx);
                if (p > 2) {
                    x_idx = std::abs(15 * ((n / 64) % 2) - x_idx);
                    y_idx = std::abs(15 * ((n / 128) % 2) - y_idx);
                }
            }
        }

        std::vector<int> si_values;
        for (int i = -A; i <= A; i += 2) {
            si_values.push_back(i);
        }

        x = si_values[x_idx];
        y = si_values[y_idx];
    }
};





class Channel {
private:
    RandomGenerator rng;

public:
    void add_noise(double& x, double& y, double gamma) {
        double noise = std::sqrt(1 / (2 * gamma));
        x += rng.rand_normal() * noise;
        y += rng.rand_normal() * noise;
    }
};


class QAMDemodulator {
private:
    int p;
    int A;

public:
    QAMDemodulator(int p_val) : p(p_val) {
        A = 2 * std::pow(2, p) - 1;
    }

    int demodulate(double x, double y) {
        int result = 0;

        x = x / 2.0 + (A / 2.0);
        y = y / 2.0 + (A / 2.0);

        x = std::round(x);
        y = std::round(y);

        x = std::max(0.0, std::min(static_cast<double>(A), x));
        y = std::max(0.0, std::min(static_cast<double>(A), y));

        if ((y > 0 && y < 3) || (y > 4 && y < 7) || (y > 8 && y < 11) || (y > 12 && y < 15)) {
            result |= 0x02;
        }
        if ((x > 0 && x < 3) || (x > 4 && x < 7) || (x > 8 && x < 11) || (x > 12 && x < 15)) {
            result |= 0x01;
        }

        if (p > 0) {
            if ((y > 1 && y < 6) || (y > 9 && y < 14)) {
                result |= 0x08;
            }
            if ((x > 1 && x < 6) || (x > 9 && x < 14)) {
                result |= 0x04;
            }

            if (p > 1) {
                if (y > 3 && y < 12) {
                    result |= 0x20;
                }
                if (x > 3 && x < 12) {
                    result |= 0x10;
                }

                if (p > 2) {
                    if (y > 7) {
                        result |= 0x80;
                    }
                    if (x > 7) {
                        result |= 0x40;
                    }
                }
            }
        }

        return result;
    }
};

void plot_results(const std::vector<double>& D_db, const std::vector<double>& pe, int p) {
    std::ofstream datafile("qam_results.dat");
    for (size_t i = 0; i < D_db.size(); ++i) {
        datafile << D_db[i] << " " << pe[i] << "\n";
    }
    datafile.close();

    std::ofstream plotfile("plot_script.gp");
    plotfile << "set terminal pngcairo enhanced font 'arial,10' fontscale 1.0 size 700, 500\n";
    plotfile << "set output 'DispPError.png'\n";
    plotfile << "set title 'D vs Bit Error Probability for QAM (p=" << p << ")'\n";
    plotfile << "set xlabel 'D (dB)'\n";
    plotfile << "set ylabel 'Bit Error Probability'\n";
    plotfile << "set logscale y\n";
    plotfile << "set grid\n";
    plotfile << "set key top right\n";
    plotfile << "plot 'qam_results.dat' with linespoints title 'Simulated Bit Error', \\\n";

    if (p == 0) {
        
        plotfile << "     0.5*erfc(sqrt(0.5/(10**(x/10)))) title 'Theoretical BER (p=0)'\n";
    }
    else {
        plotfile << "     " << (std::pow(2, p + 1) - 1) << "*erfc(sqrt(0.5/(10**(x/10))))*(2**(" << p << "+1)-" << (std::pow(2, p + 1) - 1) << "*erfc(sqrt(0.5/(10**(x/10))))/2)/(4**" << (p + 1) << ") title 'Theoretical Upper Bound', \\\n";
        plotfile << "     " << (std::pow(2, p + 1) - 1) << "*erfc(sqrt(0.5/(10**(x/10))))*(2**(" << p << "+1)-" << (std::pow(2, p + 1) - 1) << "*erfc(sqrt(0.5/(10**(x/10))))/2)/(4**" << (p + 1) << ")/" << (p + 1) << " title 'Theoretical Lower Bound', \\\n";
    }
    plotfile.close();
}

void simulate_qam(int p, int test_points, const std::vector<double>& D_db) {
    QAMModulator_simple modulator(p);
    Channel channel;
    QAMDemodulator demodulator(p);

    std::vector<double> pe(D_db.size(), 0.0);
    RandomGenerator rng;

    for (size_t j = 0; j < D_db.size(); ++j) {
        double D = db_to_linear(D_db[j]);
        double gamma = 1 / (2 * D);
        int error_count = 0;

        for (int i = 0; i < test_points; ++i) {
            int symbol = rng.rand_int(0, modulator.get_q() - 1);

            double x, y;
            modulator.modulate(symbol, x, y);
            channel.add_noise(x, y, gamma);

            int received = demodulator.demodulate(x, y);

            std::vector<int> tx_bits = int_to_bits(symbol, (p + 1) * 2);
            std::vector<int> rx_bits = int_to_bits(received, (p + 1) * 2);

            for (int k = 0; k < (p + 1) * 2; ++k) {
                if (tx_bits[k] != rx_bits[k]) {
                    error_count++;
                }
            }
        }
        pe[j] = static_cast<double>(error_count) / (test_points * (p + 1) * 2);
    }
    plot_results(D_db, pe, p);
}

int main() {
    int p = 3;
    int Test = 10000;

    std::vector<double> D_db;
    for (double g = -15.0; g <= 15.0; g += 0.2) {
        D_db.push_back(g);
    }

    simulate_qam(p, Test, D_db);
    std::cout << "plot_script created";
    return 0;
}