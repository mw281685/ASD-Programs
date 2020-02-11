#include<cstdio>
#include<cstdlib>
#include<algorithm>
#include <iostream>
#include <vector>

using namespace std;
#define DEBUG 0

int n, m;
char initialSetting[1000001];
int mod = 1000000007;

struct Macierz{
	unsigned long long T[2][2];
	Macierz(){
		T[0][0] = 0;
		T[0][1] = 0;
		T[1][0] = 0;
		T[1][1] = 0;

	}

	void operator=(const Macierz &M) {
		for (int i = 0; i < 2; i++)
			for (int j =0; j < 2; j++)
				this->T[i][j]=M.T[i][j];
	}

	void wypisz() {
		printf("%llu %llu\n%llu %llu\n", T[0][0], T[0][1], T[1][0], T[1][1]);
	}
};

struct Node{
	Node *l; // reference to left child
	Node *r; //reference to right child
	Node *p; // reference to parent
	int v; //wartosc indeksu
	int size;
	int ifRotated;
	Macierz *tablica;  // R- 0
	int RR = 0;
	int RG = 0;
	int GR = 0;
	int GG = 0;
	Node(){};
	Node(Node* _l,Node* _r, Node* _p, int _v, int _size, int _ifRotated, int _RR, int _RG, int _GR, int _GG) : l(_l), r(_r), p(_p), v(_v), size(_size), ifRotated(_ifRotated), RR(_RR), RG(_RG), GR(_GR), GG(_GG) {};

	void wypisz() {
			printf("ja=%15lu l=%15lu r=%15lu p=%15lu v=%d size=%d ifRotated=%d, RR = %d, RG = %d , GR =%d, GG = %d \n",
			(size_t)this, (size_t)l, (size_t)r, (size_t)p, v, size, ifRotated, tablica->T[0][0], tablica->T[0][1],tablica->T[1][0], tablica->T[1][1]);
			if (l != NULL) l->wypisz();
			if (r != NULL) r->wypisz();
	}
};





void Inorder(Node *R) // wypisz drzewo w inorder.
{
	if(!R) return;
	Inorder(R->l);
	if(R->p)printf("v: %d, size: %d, ifRotated: %d, p = %d   ",R->v, R->size, R->ifRotated, R->p->v);
	else printf("v: %d, size: %d, ifRotated: %d, ROOT  ",R->v, R->size, R->ifRotated);
	if(R->l) printf("l: %d,  ",R->l->v );
	if(R->r) printf("r: %d ",R->r->v);
	puts("");
	Inorder(R->r);
}

void printSequence_(Node *R, bool rotated){
	if(!R) return;
	if(R->ifRotated) rotated = !rotated;

	if (rotated) {
		printSequence_(R->r, rotated);
		printf("%c", initialSetting[R->v] );
		printSequence_(R->l, rotated);
	} else {
		printSequence_(R->l, rotated);
		printf("%c", initialSetting[R->v] );
		printSequence_(R->r, rotated);
	}
}

void printSequence(Node *R){
	printSequence_(R, false);
	puts("");
}




Macierz polaczMacierze(Macierz *A, Macierz  *B){
	Macierz wynik;
	wynik.T[0][0] = (A->T[0][0]+
					   B->T[0][0]+
					 A->T[0][0]*(B->T[1][0] + B->T[0][0])%mod+
					 A->T[0][1]*(B->T[0][0])%mod)%mod;

	wynik.T[0][1] = (A->T[0][1] +
					B->T[0][1] +
			(A->T[0][0]*(B->T[1][1] + B->T[0][1]))%mod +
			(A->T[0][1]*B->T[0][1]))%mod;


	wynik.T[1][0] = (A->T[1][0] +
					  B->T[1][0] +
			(A->T[1][0]*(B->T[0][0] + B->T[1][0]) + A->T[1][1]*B->T[0][0])%mod)%mod; // MS modified: + A->T[1][1]*B->T[0][0]

	wynik.T[1][1] = (A->T[1][1] +
			B->T[1][1]+
			(A->T[1][0]*(B->T[0][1] + B->T[1][1]))%mod+
			(A->T[1][1]*(B->T[0][1]))%mod)%mod;
	return wynik;
}

void przeliczWyn(Node* P){
	if (P==NULL) return;
	//clear(P);
	//printf("Jestem w fun PrzeliczWyn \n");
	Macierz lM, sM, rM;

	// inicjalizacja noda
	if(initialSetting[P->v] == 'R' ){
		sM.T[0][0] = 1;
		sM.T[0][1] = 0;
		sM.T[1][0] = 0;
		sM.T[1][1] = 0;
	}else{
		sM.T[0][0] = 0;
		sM.T[0][1] = 0;
		sM.T[1][0] = 0;
		sM.T[1][1] = 1;
	}

	if(P->l){ // jest lewy syn
		lM.T[0][0] = P->l->tablica->T[0][0];
		lM.T[0][1] = P->l->tablica->T[0][1];
		lM.T[1][0] = P->l->tablica->T[1][0];
		lM.T[1][1] = P->l->tablica->T[1][1];
	}

	if(P->r){ // jest prawy syn
		rM.T[0][0] = P->r->tablica->T[0][0];
		rM.T[0][1] = P->r->tablica->T[0][1];
		rM.T[1][0] = P->r->tablica->T[1][0];
		rM.T[1][1] = P->r->tablica->T[1][1];
	}

	Macierz wynL;
	if(P->l){
		//clear(P->l);
		wynL = polaczMacierze(&lM, &sM);
	}else{
		wynL = sM;
	}

	Macierz wynR;
	if(P->r){
		//clear(P->r);
		wynR = polaczMacierze(&wynL, &rM);
	}else wynR = wynL;


// update liczbe ciagow dla noda
	P->tablica->T[0][0] = wynR.T[0][0];
	P->tablica->T[0][1] = wynR.T[0][1];
	P->tablica->T[1][0] = wynR.T[1][0];
	P->tablica->T[1][1] = wynR.T[1][1];
//	cout<< "dla P->v "<< P->v << " P RR to "<< P->RR <<endl;
}



void swapSonsIfRotated(Node *P)
{
	if(!P) return;
	Node* tmp = P->l;
	P->l = P->r;
	P->r = tmp;
	P->ifRotated = 1- (P->ifRotated);

	if(P->l){
		P->l->ifRotated = 1- (P->l->ifRotated);
		int tmpValL = P->l->tablica->T[0][1];
		P->l->tablica->T[0][1] = P->l->tablica->T[1][0];
		P->l->tablica->T[1][0] = tmpValL;
	}
	if(P->r){
		P->r->ifRotated = 1 -(P->r->ifRotated);
		int tmpValR = P->r->tablica->T[0][1];
		P->r->tablica->T[0][1] = P->r->tablica->T[1][0];
		P->r->tablica->T[1][0] = tmpValR;
	}

	// update ( zamien RG z GR!!!) :
	//int tmpVal = P->tablica->T[0][1];
	//P->tablica->T[0][1] = P->tablica->T[1][0];
	//P->tablica->T[1][0] = tmpVal;

	przeliczWyn(P); // MS added 11/02/2020!!!
}


void clear(Node *P){
		if(P == NULL) return;
		if(P->ifRotated){
			swapSonsIfRotated(P);
			if(P->l && P->l->ifRotated) swapSonsIfRotated(P->l);
			if(P->r && P->r->ifRotated) swapSonsIfRotated(P->r);
		}
}






void updateSize(Node* P)
{
	//if (P!=NULL) printf("Wchodze  do funkcji przelicz Wyn dla node %d \n", P->v);
	clear(P);
	if(P->l && P->r){
		P->size = P->l->size + P->r->size + 1;
	}
	else if (P->l){
		P->size  = P->l->size + 1;
	}
	else if (P->r ){
		P->size = P->r->size + 1;
	}
	else P->size = 1;
	P = P->p;
	clear(P);

	przeliczWyn(P);
	//if (P!=NULL) printf("Wyszlam z funkcji przelicz Wyn dla node %d \n", P->v);

}


void rightRotate(Node *P)
{
//	if(P->ifRotated) swapSonsIfRotated(P); // ms commented
	clear(P);
	Node *T=P->l;
	Node *B=T->r;
	Node *D=P->p;
	if(D)
	{
		clear(D);
		if(D->r==P) D->r=T;
		else D->l=T;
	}
	if(B){
		clear(B);
		B->p=P;
	}
	T->p=D;
	T->r=P;

	P->p=T;
	P->l=B;
	updateSize(P);
	przeliczWyn(P);
	updateSize(T);
	przeliczWyn(T);
}



void leftRotate(Node *P)
{
	//if(P->ifRotated) swapSonsIfRotated(P);
	clear(P);
	Node *T=P->r;
	Node *B=T->l;
	Node *D=P->p;
	if(D)
	{
		clear(D);
		if(D->r==P) D->r=T;
		else D->l=T;
	}
	if(B){
		clear(B);
		B->p=P;
	}

	T->p=D;
	T->l=P;

	P->p=T;
	P->r=B;
	updateSize(P);
	przeliczWyn(P);
	updateSize(T);
	przeliczWyn(T);
}

Node* Splay(Node *T){
	//if(T->p && T->p->p) printf("Robie splay na %d , p : %d , pp = %d \n", T->v, T->p->v, (T->p)->p->v);
	clear(T);
	while(true)
	{
		Node *p=T->p; //PARENT
		if(!p) break;
		Node *pp=p->p; // GRANDPARENT
		if(!pp)//Zig
		{
			if(p->l==T) rightRotate(p); // && (T->ifRotated == 0 ))||(p->r == T && (T->ifRotated == 1)
			else leftRotate(p);
			break;
		}
		if(pp->l==p)
		{
			if(p->l==T)
			{//ZigZig
				rightRotate(pp);
				rightRotate(p);
			}
			else
			{//ZigZag
				leftRotate(p);
				rightRotate(pp);
			}
		}
		else
		{
			if(p->l==T)
			{//ZigZag
				rightRotate(p);
				leftRotate(pp);
			}
			else
			{//ZigZig
				leftRotate(pp);
				leftRotate(p);
			}
		}
	}
	return T;
}


Node* Insert(int n){

	Node* node =(Node *)malloc(sizeof(Node));
	node->tablica = new Macierz();

	node->v = 1;
	node->ifRotated = 0;
	node->l = NULL;
	node->r = NULL;
	node-> p = NULL;
	node->size = 1;

	if(initialSetting[1] == 'G'){
		node->tablica->T[1][1]=1;
	}else{
		node->tablica->T[0][0]=1;
	}

	for(int i = 2; i <= n ; i++){
		Node* n2= (Node *)malloc(sizeof(Node));
		Macierz *tab =(Macierz*)malloc(sizeof(Macierz));
		n2->tablica=tab;
		node->p = n2;
		n2->v = i;
		n2->ifRotated = 0;
		n2->l = node;
		n2->r = NULL;
		n2->p = NULL;
		n2->size = 1 + node->size;
		if(initialSetting[i] == 'G'){ // mamy samych lewych synow;
			n2->tablica->T[1][1]=n2->l->tablica->T[1][0]+n2->l->tablica->T[1][1]+1;
			n2->tablica->T[1][0]=n2->l->tablica->T[1][0];
			n2->tablica->T[0][1]=n2->l->tablica->T[0][0]+n2->l->tablica->T[0][1];
			n2->tablica->T[0][0]=n2->l->tablica->T[0][0];

		}else{
			n2->tablica->T[1][1]=n2->l->tablica->T[1][1];
			n2->tablica->T[1][0]=n2->l->tablica->T[1][1]+n2->l->tablica->T[1][0]*2;
			n2->tablica->T[0][1]=n2->l->tablica->T[0][1];
			n2->tablica->T[0][0]=n2->l->tablica->T[0][1]+n2->l->tablica->T[0][0]*2+1;
		}
		node = n2;
	}
	return node;
}


Node* Find(Node* root, int idx) // daje wskaznik na wierzcholek w kolejnosci idx, lub NULL jak takiego nie ma
{
	clear(root);
	int idxOriginal = idx;
	if(!root) return NULL;
	Node *P=root;

	int sizeL =0, sizeR = 0;
	while(P)
	{
		//if(P->ifRotated) swapSonsIfRotated(P);
		clear(P);
		sizeR = (P->r)? P->r->size : 0;
		sizeL = (P->l)? P->l->size : 0;

		if (sizeL + 1 == idx){
			break;
		}

		if(idx<= sizeL)
		{
			if(P->l){
				clear(P->l);
				P=P->l;
			}
			else
				break;
		}
		else
		{
			if(P->r){
				clear(P->r);
				P=P->r;
				idx = idx - sizeL - 1;
			}
			else
				break;
		}
	}

	root = Splay(P);

	return P;

}


pair<Node*, Node*> split(Node* root, int idx, bool ifLeft)
{

	if(DEBUG){
		printf("Jestem w split, root->v = %d , idx = %d \n", root->v, idx);
		printf("Drzewo przed operacja find \n");
		root->wypisz();
	}

	Node* n = Find(root, idx); // robi splay na node = idx
	clear(n);

	if(DEBUG){
		printf("Root po odnalezieniu wierzcholka do splitowania : root->v = %d :\n ", n->v);
	  n->wypisz();
		Inorder(n);
		printf("-------------\n");
	}


	if(n && ifLeft && n->l ){
		Node* leftSplit = n->l;
		clear(leftSplit);
		n->size -=n->l->size;
		n->l = NULL;
		leftSplit->p = NULL;
		przeliczWyn(n);
		return make_pair(n, leftSplit);
	}else if(n && !ifLeft && n->r){
		Node* rightSplit = n->r;
		clear(rightSplit);
		rightSplit->p = NULL;
		n->size -=n->r->size;
		n->r = NULL;
		przeliczWyn(n);
		return make_pair(n, rightSplit);
	}
	else return make_pair(n, (Node*)NULL);

}


	Node* merge(Node* T1, Node* T2){
	if(!T1 && !T2 ) return NULL;
	else if (!T1){
		clear(T2);
		return T2;
	}
	else if(!T2){
		clear(T1);
		return T1;
	}

	if(T1->ifRotated) swapSonsIfRotated(T1);
	while(T1->r){
		T1 = T1->r;
		if(T1->ifRotated){
			swapSonsIfRotated(T1);
		}
	}
	if(T1->ifRotated) swapSonsIfRotated(T1);
	Splay(T1);
	T2->p = T1;
	T1->r = T2;
	updateSize(T1);
	przeliczWyn(T1);
	return T1;
}

Node* reverseSplay(Node* root, int leftIdx, int rightIdx){

	if(DEBUG){
		printf("-----------Reverse : Wyjsciowe drzewo ------------ \n");
		root->wypisz();
		printf("------------------\n");
		printf("Reverse leftIdx = %d , rightIdx = %d \n", leftIdx, rightIdx);
	}
	clear(root);
	pair<Node*, Node*> x;
	x = split(root, leftIdx, 1);
	Node* leftTree = x.second;
	clear(leftTree);
	root = x.first;
	clear(root);

	if(DEBUG){
		printf("Wypisz lewe drzewo \n");
		if (leftTree) leftTree->wypisz();
		printf("------\n");

		printf("root po oddzieleniu lewego drzewa: \n");
		if (root) {
			root->wypisz();
			Inorder(root);
		}
		printf("------\n");

	}

   	x = split(root, rightIdx - leftIdx  + 1 , 0);

	root = x.first;
	clear(root);
	Node* rightTree = x.second;
	clear(rightTree);

	if(DEBUG){
		printf("Wypisz prawe drzewo  \n");
		if (rightTree) Inorder(rightTree);
		printf("------\n");
	}

	root->ifRotated = !root->ifRotated;
	clear(root);
	if(DEBUG){
		printf("Wypisz srodkowe zrotowane drzewo  \n");
		if (root) root->wypisz();
		printf("------\n");
	}


	Node* mergeLeft = merge(leftTree, root);
	clear(mergeLeft);

	//printf("Wypisz mergeLeft  \n");
	//Inorder(mergeLeft);
	//printf("------\n");

	//printf("Wypisz rightTree\n");
	//Inorder(rightTree);
	//printf("------\n");

	Node* mergeRight = merge(mergeLeft, rightTree);
	clear(mergeRight);
/*
	printf("Wypisz mergeRight  \n");
	Inorder(mergeRight);
	printf("------\n");
*/
	root = mergeRight;

	if(DEBUG){
		printf("Wypisz drzewo po reverse!  \n");
		root->wypisz();
		printf("------------\n");
	}
/*	Inorder(root);
	printf("------\n");
*/
	return root;
}

pair<Node*, long long> calculateResult(Node* root, int leftIdx, int rightIdx){
	clear(root);
	pair<Node*, Node*> x;
	x = split(root, leftIdx, 1);
	Node* leftTree = x.second; // drzewo do odrzucenia
	clear(leftTree);
	root = x.first;
	clear(root);

	if(DEBUG){
		printf("Po splay na leftIdx = %d drzewo srodkowe to : \n", leftIdx);
		root->wypisz();
		printf("---------------- \n");
	}

  x = split(root, rightIdx - leftIdx  + 1 , 0);
	root = x.first;
	clear(root);

	if(DEBUG){
		printf("Po splay na rightIdx = %d drzewo srodkowe to : \n", rightIdx);
		root->wypisz();
		printf("---------------- \n");
	}


	Node* rightTree = x.second;
	clear(rightTree);
	long long wyn = (root->tablica->T[0][0]+root->tablica->T[0][1]+root->tablica->T[1][0]+root->tablica->T[1][1])%mod;

	Node* nowy = merge(leftTree, root);
	clear(nowy);
	Node* nowszy = merge(nowy, rightTree);
	clear(nowszy);
	return make_pair(nowszy, wyn);
	//return (((root->RR + root->RG)%mod + root-> GR)%mod + root-> GG)%mod;
}


int main(){
	char buf[20];
	scanf("%d%d", &n, &m);
	scanf("%s", initialSetting+1);

  //for(int i=1; i<=n; i++) printf("Armia[%d] = %c \n", i, initialSetting[i]);
	Node *root = Insert(n);

	if(DEBUG){
		printf("Wypisuje drzewo po inicjacji: : \n");
		root->wypisz();
		printf("--------\n");
	}

	int a, b;

	for (int i = 0; i < m; ++i) {
		scanf("%s%d%d", buf, &a, &b);
		if (buf[0] == 'O') {
			root = reverseSplay(root, a, b);
			clear(root);
			//printf("Ostateczna kolejnosc : ---- \n");
//			printSequence(root);
			//Inorder(root);
		} else{
			pair<Node*, long long> wynik = calculateResult(root, a, b);
			//printf("wynik to:");
			printf("%lld\n", (wynik.second)%mod);

			root = wynik.first;
			clear(root);
			if(DEBUG){
				printf("Drzewo po wykonaniu zapytania \n");
				root->wypisz();
				printf("$$$$$$$$$$$$$$$$$$$$\n");
			}
		}
	}

	//printSequence(root);
	//printf("\n");
}
