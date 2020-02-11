#include<cstdio>
#include<cstdlib>
#include<algorithm>

using namespace std;

int main(int argc, char** argv){
  srand(atoi(argv[1]));
  int n = 5 + rand()%3;
  int m = 5 + rand()%3;

  printf("%d %d\n",n, m);

  for(int i = 0; i< n; i++ ){
    if (rand()%100 < 50) printf("R");
    else printf("G");
  }
  printf("\n");

  for(int i=0; i<m; i++){
    int a = rand()%n + 1;
    int b = rand()%n + 1;
    if(a > b) swap(a, b);

    if(rand()%100 < 50){
        printf("O %d %d\n", a, b );
    }else printf("? %d %d \n", a, b);
  }

  return 0 ;
}
