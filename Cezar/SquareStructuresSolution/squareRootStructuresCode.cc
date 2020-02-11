#include <cstdio>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;

const int P = 1000000007;

struct macierz {
	unsigned int t[3][3];
	macierz() {}

	void wyczysc() {
		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j) {
				t[i][j] = 0;
			}
		}
	}
	void wypisz() {
		for (int i = 0; i < 3; i++) {
			for (int j =0; j < 3; ++j) {
				printf("%d ", t[i][j]);
			}
			printf("\n");
		}
	}

	bool operator==(const macierz &M) const {
		for (int i = 0; i < 3; i++) for (int j =0; j < 3; j++) if (M.t[i][j] != t[i][j]) return false;
		return true;
	}
};

struct wektor {
	unsigned int t[3];
	wektor() { t[0] = t[1] = t[2] = 0; }
};

inline unsigned int modula(long long x) {
	if (x < P) return x; else return x % P;
}

macierz iloczyn(const macierz &A, const macierz &B) {
	macierz wynik;
	for (int i = 0; i < 3; ++i) {
		wynik.t[i][0] = modula((long long)A.t[i][0] * B.t[0][0] +(long long)A.t[i][1] * B.t[1][0] + (long long)A.t[i][2] * B.t[2][0]);
		wynik.t[i][1] = modula((long long)A.t[i][0] * B.t[0][1] +(long long)A.t[i][1] * B.t[1][1] + (long long)A.t[i][2] * B.t[2][1]);
		wynik.t[i][2] = modula((long long)A.t[i][0] * B.t[0][2] +(long long)A.t[i][1] * B.t[1][2] + (long long)A.t[i][2] * B.t[2][2]);
	}
	return wynik;
}

wektor iloczyn(const wektor &w, const macierz &B) {
	wektor wynik;
	wynik.t[0] = modula((long long)w.t[0] * B.t[0][0] +(long long)w.t[1] * B.t[1][0] + (long long)w.t[2] * B.t[2][0]);
	wynik.t[1] = modula((long long)w.t[0] * B.t[0][1] +(long long)w.t[1] * B.t[1][1] + (long long)w.t[2] * B.t[2][1]);
	wynik.t[2] = modula((long long)w.t[0] * B.t[0][2] +(long long)w.t[1] * B.t[1][2] + (long long)w.t[2] * B.t[2][2]);
	
	return wynik;
}
struct przedzial {
	int p,k;
	bool odw;
	przedzial() {}
	przedzial(int _p, int _k, bool _odw) : p(_p), k(_k), odw(_odw) {}
};

int n, m;
const int max_n = 100100;
char z[max_n], z2[max_n];
vector<przedzial> v;
macierz gal, rzymianin, gal_odw, rzymianin_odw;
macierz ilo_lewy[max_n], ilo_lewy_odw[max_n], ilo_prawy[max_n], ilo_prawy_odw[max_n];
macierz orig_ilo_lewy[max_n], orig_ilo_lewy_odw[max_n]; 
macierz krok[1<<8], krok_odw[1<<8];
bool lewy_init[max_n], lewy_odw_init[max_n], prawy_init[max_n], prawy_odw_init[max_n];

void inicjalizuj_v() {
	v.clear();
	v.push_back(przedzial(1, n, false));
	macierz jednostkowa;
	jednostkowa.wyczysc();
	jednostkowa.t[0][0] = jednostkowa.t[1][1] = jednostkowa.t[2][2] = 1;
	ilo_lewy[0] = ilo_lewy_odw[0] = jednostkowa;
	orig_ilo_lewy[0] = orig_ilo_lewy_odw[0] = jednostkowa;
	for (int i = 1; i <= n; ++i) {
		lewy_init[i] = lewy_odw_init[i] = false;
		prawy_init[i] = prawy_odw_init[i] = false;
		//orig_ilo_lewy[i] = iloczyn(orig_ilo_lewy[i-1], (z[i] == 'R') ? rzymianin : gal);
		//orig_ilo_lewy_odw[i] = iloczyn((z[i] == 'R') ? rzymianin_odw : gal_odw, orig_ilo_lewy_odw[i-1]);
	}
	lewy_init[0] = lewy_odw_init[0] = true;

	for (int i = 1; i <= n / 8; ++i) {
		int skok = 0;
		for (int j = 0; j < 8; ++j) skok = skok * 2 + ((z[(i-1)*8+j+1] == 'R') ? 1 : 0);
		//printf("i=%d skok=%d\n", i,skok);
		ilo_lewy[i * 8] = iloczyn(ilo_lewy[(i-1)*8], krok[skok]);
		ilo_lewy_odw[i * 8] = iloczyn(krok_odw[skok], ilo_lewy_odw[(i-1)*8]);
		lewy_init[i*8] = lewy_odw_init[i*8] = true;
	}

	ilo_prawy[n+1] = ilo_prawy_odw[n+1] = jednostkowa;
	//fprintf(stderr, "init prawy\n");
	for (int i = 1; i <= n / 8; ++i) {
		int skok = 0;
		for (int j = 0; j < 8; ++j) skok = skok * 2 + ((z[n-(i-1)*8-j] == 'R') ? 1 : 0);	
		//fprintf(stderr, "i=%d skok=%d\n", i, skok);
		ilo_prawy[n+1-i*8] = iloczyn(ilo_prawy[n+1 - (i-1) *8], krok[skok]);
		ilo_prawy_odw[n+1-i*8] = iloczyn(krok_odw[skok], ilo_prawy_odw[n+1-(i-1)*8]);
		prawy_odw_init[n+1-i*8] = prawy_init[n+1-i*8] = true;
	}
	//fprintf(stderr, "init prawy koniec\n");
	
	/*for (int i=n; i >= 1; --i) {
		ilo_prawy[i] = iloczyn(ilo_prawy[i+1], (z[i] == 'R')? rzymianin : gal);
		ilo_prawy_odw[i] = iloczyn((z[i] == 'R') ? rzymianin_odw : gal_odw, ilo_prawy_odw[i+1]);
	}*/
}

void inicjalizuj_macierze() {
	rzymianin.wyczysc();
	gal.wyczysc();
	rzymianin.t[0][0] = 2;
	rzymianin.t[1][0] = rzymianin.t[1][1] = 1;
	rzymianin.t[2][0] = rzymianin.t[2][2] = 1;
	gal.t[0][0] = gal.t[0][1] = 1;
	gal.t[1][1] = 1;
	gal.t[2][1] = gal.t[2][2] = 1;

	rzymianin_odw.wyczysc();
	gal_odw.wyczysc();
	rzymianin_odw.t[0][0] = (P + 1) / 2;
	rzymianin_odw.t[1][0] = rzymianin_odw.t[2][0] = (P - 1) / 2;
	rzymianin_odw.t[1][1] = rzymianin_odw.t[2][2] = 1;
	gal_odw.t[0][0] = gal_odw.t[1][1] = gal_odw.t[2][2] = 1;
	gal_odw.t[0][1] = gal_odw.t[2][1] = P - 1;

	macierz jednostkowa;
	jednostkowa.t[0][0] = jednostkowa.t[1][1] = jednostkowa.t[2][2] = 1;

	for (int i = 0; i < (1<<8); ++i) {
		krok[i] = krok_odw[i] = jednostkowa;
		for (int j = 7; j >= 0; --j) {
			krok[i] = (i&(1<<j)) ? iloczyn(krok[i], rzymianin) : iloczyn(krok[i], gal);
			krok_odw[i] = (i&(1<<j)) ? iloczyn(rzymianin_odw, krok_odw[i]) : iloczyn(gal_odw, krok_odw[i]);		
		}
	}
}

void porzadek() {
	int index = 0;
	for (int i = 0; i < v.size(); ++i) {
		if (!v[i].odw) {
			for (int j = v[i].p; j <= v[i].k; ++j) z2[++index] = z[j];
		} else {
			for (int j = v[i].k; j >= v[i].p; --j) z2[++index] = z[j];
		}
	}
	for (int i = 1; i <= n; ++i) z[i] = z2[i];
	inicjalizuj_v();
}

void odwroc(int pocz, int kon) {
	int suma = 0;
	int bez_zmian = 0;
	vector<przedzial> lewy, srodek, prawy;
	for (int i = 0; i < v.size(); ++i) {
		int p = pocz - suma;
		int k = kon - suma;
		int d = v[i].k - v[i].p + 1;
	
		if (p > d) {
			++bez_zmian;
		} else if (k < 1) {
			int old_v_size = v.size();
			v.resize(bez_zmian + lewy.size() + srodek.size() + prawy.size() + v.size() - i);
			for (int j = old_v_size - 1; j >= i; --j) v[j + (v.size() - old_v_size)] = v[j];
			copy(lewy.begin(), lewy.end(), v.begin() + bez_zmian);
			int pozycja = lewy.size() + bez_zmian;
			reverse(srodek.begin(), srodek.end());
			copy(srodek.begin(), srodek.end(), v.begin() + pozycja);
			pozycja += srodek.size();
			copy(prawy.begin(), prawy.end(), v.begin() + pozycja);
			return;
			//prawy.emplace_back(v[i]);
		} else {
			if (!v[i].odw) {
				if (p > 1) {
					lewy.emplace_back(przedzial(v[i].p, v[i].p + p - 2, v[i].odw));
				}
				srodek.emplace_back(przedzial(max(v[i].p + p - 1, v[i].p), min(v[i].k, v[i].p + k - 1), !v[i].odw));
				if (v[i].p + k - 1 < v[i].k) {
					prawy.emplace_back(przedzial(v[i].p + k, v[i].k, v[i].odw));
				}
			} else {
				if (p > 1) {
					lewy.emplace_back(przedzial(v[i].k - (p - 2), v[i].k, v[i].odw));
				}
				int a = max(1, p);
				int b = min(d, k);
				srodek.emplace_back(przedzial(v[i].k - b + 1, v[i].k - a + 1, !v[i].odw));
				if (d > k) {
					prawy.emplace_back(przedzial(v[i].p, v[i].p + d - k - 1, v[i].odw));
				}
			}
		}
		suma += v[i].k - v[i].p + 1;
	}
	v.resize(bez_zmian + lewy.size() + srodek.size() + prawy.size());
	copy(lewy.begin(), lewy.end(), v.begin() + bez_zmian);
	int pozycja = lewy.size() + bez_zmian;
	reverse(srodek.begin(), srodek.end());
	copy(srodek.begin(), srodek.end(), v.begin() + pozycja);
	pozycja += srodek.size();
	copy(prawy.begin(), prawy.end(), v.begin() + pozycja);
}

macierz wynik_przedzial(int p, int k, bool odw) {
	//fprintf(stderr, "wynik_przedzial %d %d %d\n", p, k, odw);
	if (!odw) {
		if (!lewy_odw_init[p-1]) {
			int a = ((p-1)/8) * 8;
			while (a < p-1) {
				++a;
				ilo_lewy_odw[a] = iloczyn((z[a] == 'R') ? rzymianin_odw : gal_odw, ilo_lewy_odw[a-1]);
				lewy_odw_init[a] = true;
			}
		}
		if (!lewy_init[k]) {
			int a = (k/8) * 8;
			while (a < k) {
				++a;
				ilo_lewy[a] = iloczyn(ilo_lewy[a-1], (z[a] == 'R') ? rzymianin : gal);
				lewy_init[a] = true;
			}
		}

		return iloczyn(ilo_lewy_odw[p-1], ilo_lewy[k]);
	} else {
		if (!prawy_odw_init[k+1]) {
			int a = n+1 - ((n+1 - (k+1))/8) * 8;
			while (a > k + 1) {
				--a;
				ilo_prawy_odw[a] = iloczyn((z[a] == 'R') ? rzymianin_odw : gal_odw, ilo_prawy_odw[a+1]);
				prawy_odw_init[a] = true; 
			}
		}
		if (!prawy_init[p]) {
			int a = n+1 - ((n+1 - p)/8) * 8;
			while (a > p) {
				--a;
				ilo_prawy[a] = iloczyn(ilo_prawy[a+1], (z[a] == 'R')? rzymianin : gal);
				prawy_init[a] = true;
			}
		}

		return iloczyn(ilo_prawy_odw[k+1], ilo_prawy[p]);
	}
}

int policz_wynik(int pocz, int kon) {
	int suma = 0;
	macierz jednostkowa;
	wektor wynik;
	jednostkowa.wyczysc();
	jednostkowa.t[0][0] = jednostkowa.t[1][1] = jednostkowa.t[2][2] = 1;
	wynik.t[0] = wynik.t[1] = 0;
	wynik.t[2] = 1;

	for (int i = 0; i < v.size(); ++i) {
		int p = pocz - suma;
		int k = kon - suma;
		if (k <= 0) break;
		int d = v[i].k - v[i].p + 1;
		//if (d < p) continue;
		if (p <= 1 && d <= k) {
		//	wynik_przedzial(v[i].p, v[i].k, v[i].odw).wypisz();
			wynik = iloczyn(wynik, wynik_przedzial(v[i].p, v[i].k, v[i].odw));
		} else {
			int a = max(1, p);
			int b = min(k, d);	
			if (a <= b) {
				if (!v[i].odw) {
					wynik = iloczyn(wynik, wynik_przedzial(v[i].p + a - 1, v[i].p + b - 1, v[i].odw));
				} else {
					wynik = iloczyn(wynik, wynik_przedzial(v[i].k - b + 1, v[i].k - a + 1, v[i].odw));
				}
			}
		}
		suma += d;
	}
	//wynik.wypisz();
	return (wynik.t[0] + wynik.t[1]) % P;	
}

void wypisz() {
	//for (int i = 0; i < v.size(); ++i) {
	//	printf("%d %d %d\n", v[i].p, v[i].k, v[i].odw);
	//}
	for (int i = 0; i < v.size(); ++i) {
		if (!v[i].odw) {
			for (int j = v[i].p; j <= v[i].k; ++j) printf("%c", z[j]);
		} else {
			for (int j = v[i].k; j >= v[i].p; --j) printf("%c", z[j]);
		}
	}
	printf("\n");
}

int main() {
	int a,b;
	char buf[20];
	scanf("%d%d", &n, &m);
	scanf("%s", z+1);
	inicjalizuj_macierze();
	inicjalizuj_v();
/*	for (int i = 1; i <= n; i++) {
		printf("lewy=%d\n", i);
		iloczyn(ilo_lewy[i], ilo_lewy_odw[i]).wypisz();
	}
	for (int i = 1; i <= n; i++) {
		printf("prawy=%d\n", i);
		iloczyn(ilo_prawy[i], ilo_prawy_odw[i]).wypisz();
	}*/
	for (int i = 0; i < m; ++i) {
		scanf("%s%d%d", buf, &a, &b);
		if (buf[0] == 'O') {
			odwroc(a, b);
			if (v.size() > 2 * sqrt(n)) porzadek();
		} else {
			printf("%d\n", policz_wynik(a, b));
		}
		//wypisz();
	}
return 0;
}