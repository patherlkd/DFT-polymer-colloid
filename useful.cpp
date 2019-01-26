#include "useful.h"
#include <iostream>

using namespace std;

db Htheta(db x) {

    if (x < 0.0)
        return 0.0;
    else
        return 1.0;
}

db heaviside(db R, db x) {
    x = abs(x);
    if (x > R) {
        return 0;
    } else {
        return 1.0;
    }
}



db pow12(db g) {

    return g * g * g * g * g * g * g * g * g * g * g*g;
}



void modulate(vec &v, int rows) {

    for (int k = 0; k < rows; k++) {
        if (v(k) < 0) {
            v(k) = -v(k);
        }
    }
}

void modulate1(vec &v, int rows) {

    for (int k = 0; k < rows; k++) {
        if (v(k) < 0) {
            v(k) = 0;
        }
    }
}

void Zero_vec(vec &v, int rows) {
    for (int i = 0; i < rows; i++) {
        v(i) = 0;
    }
}

void Zero_mat(mat &m, int rows, int cols) {

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            m(i, j) = 0;
        }
    }

}
// ========================== FOURIER CODES ==============================

void fourier(vec &data, int isign) {

    int n, mmax, m, j, istep, i;
    db wtemp, wr, wpr, wpi, wi, theta, tempr, tempi;

    int nn = (data.rows() * data.cols()) / 2.0;

    n = nn << 1;
    j = 1;
    for (i = 1; i < n; i += 2) {
        if (j > i) {
            SWAP(data(j - 1), data(i - 1));
            SWAP(data(j), data(i));
        }
        m = nn;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    //Daniel - Lanczos section begins here
    mmax = 2;
    while (n > mmax) {
        istep = mmax << 1;
        theta = isign * ((2.0 * pi) / mmax);
        wtemp = sin(0.5 * theta);
        wpr = -2.0 * sqr_d(wtemp);
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for (m = 1; m < mmax; m += 2) {
            for (i = m; i <= n; i += istep) {
                j = i + mmax;
                tempr = wr * data(j - 1)-(wi * data(j));
                tempi = wr * data(j) + (wi * data(j - 1));
                data(j - 1) = data(i - 1) - tempr;
                data(j) = data(i) - tempi;
                data(i - 1) += tempr;
                data(i) += tempi;
            }
            wr = (wtemp = wr) * wpr - wi * wpi + wr;
            wi = wi * wpr + wtemp * wpi + wi;
        }
        mmax = istep;
    }


}

void Rfourier(vec &data, int isign) {
    int i, i1, i2, i3, i4;
    db c1 = 0.5, c2, h1r, h1i, h2r, h2i, wr, wi, wpr, wpi, wtemp, theta;

    int n = data.rows() * data.cols();
    theta = pi / (db) (n >> 1);
    if (isign == 1) {
        c2 = -0.5;
        fourier(data, 1);
    } else {
        c2 = 0.5;
        theta = -theta;
    }

    wtemp = sin(0.5 * theta);
    wpr = -2.0 * wtemp*wtemp;
    wpi = sin(theta);
    wr = 1.0 + wpr;
    wi = wpi;
    for (i = 1; i < (n >> 2); i++) {
        i2 = 1 + (i1 = i + i);
        i4 = 1 + (i3 = n - i1);
        h1r = c1 * (data(i1) + data(i3));
        h1i = c1 * (data(i2) - data(i4));
        h2r = -c2 * (data(i2) + data(i4));
        h2i = c2 * (data(i1) - data(i3));
        data(i1) = h1r + wr * h2r - wi*h2i;
        data(i2) = h1i + wr * h2i - wi*h2r;
        data(i3) = h1r - wr * h2r + wi*h2i;
        data(i4) = -h1i + wr * h2i + wi*h2r;
        wr = (wtemp = wr) * wpr - wi * wpi + wr;
        wi = wi * wpr + wtemp * wpi + wi;
    }
    if (isign == 1) {
        data(0) = (h1r = data(0)) + data(1);
        data(1) = h1r - data(1);
    } else {
        data(0) = c1 * ((h1r = data(0)) + data(1));
        data(1) = c1 * (h1r - data(1));
        fourier(data, -1);
    }

}

// ===========================================================================================

void tridag(vec &a, vec &b, vec &c, vec &r, vec &sol) {
    db bet;
    int n = a.rows() * a.cols();
    vec gam(n);

    if (b(0) == 0.0) {
        cout << "ERROR_1  b(0) = 0.0;" << endl;
    }

    bet = b(0);

    sol(0) = r(0) / bet;

    for (int j = 1; j < n; j++) {
        gam(j) = c(j - 1) / bet;
        bet = b(j) - (a(j) * gam(j));

        if (bet == 0.0) {
            cout << " ERROR_2 bet == 0.0." << endl;
        }

        sol(j) = (r(j) - a(j) * sol(j - 1)) / bet;

    }

    for (int j = (n - 2); j >= 0; j--) {
        sol(j) -= gam(j + 1) * sol(j + 1);

    }

}

// CLASSES

////// ARRAY ////////////////////

dbarray::dbarray() {
    setnumelems(1);
    elems = new db[numelems];
    for (int i = 0; i < numelems; i++) {
        elems[i] = 0;
    }
}

dbarray::~dbarray() {
    delete[] elems;
}

void dbarray::divideelems(db num) {

    for (int i = 0; i < numelems; i++) {
        elems[i] /= num;
    }

}

void dbarray::setelem(int index, db val) {
    if (index > numelems) {
        cout << "ERROR: index for elems > numelems. Exiting." << endl;
        exit(1);
    } else;
    elems[index] = val;
}

void dbarray::addelem(int index) {
    if (index > numelems) {
        cout << "ERROR: index for elems > numelems. Exiting." << endl;
        exit(1);
    } else;
    elems[index]++;
}

void dbarray::setnumelems(int nm) {
    numelems = nm;
}

void dbarray::resize(int size) {
    delete[] elems;
    numelems = size;
    elems = new db[numelems];
    zero();
}

void dbarray::zero() {
    int c = 0;
    while (c < numelems) {
        elems[c] = 0;
        c++;
    }

}

void dbarray::printall() {
    int i = 0;
    while (i < numelems) {
        cout << "Elems[" << i << "] = " << elems[i] << endl;
        i++;
    }
}

db dbarray::getelem(int index) {
    return elems[index];
}
