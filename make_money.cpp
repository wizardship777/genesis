#include <bits/stdc++.h>

using namespace std;

int main() {
    /*
        APPEND: TODO/FIX by \inf+ AD:
        //int -> bigint?
    */
    for (int ns = 0; ; ns++) {
        vector <string> cur_companies_s_p_500;
        vector <string> next_companies_s_p_500;
        if (cur_companies_s_p_500 != next_companies_s_p_500) {
            continue;
        }
        vector <double> cur_prices_s_p_500, next_prices_s_p_500;
        for (int i = 0; i < 500; i++) {
            if (cur_prices_s_p_500[i] < next_prices_s_p_500[i]) {
                clone(company[i], ns);
                sell(company[i], ns + 1);
            }
        }
    }
    //WHO ABOVE?
    //APPEND above?
}
