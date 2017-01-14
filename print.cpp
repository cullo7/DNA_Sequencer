#include <iostream>

using namespace std;

void print_num(char a, int num){
  for(int i = 0;i < num;i++){
      cout << a;
  }
}

int main(int argc, char * argv[]){

  int bases[] = {'A', 'T', 'G', 'C'};
  for(int i = 0; i < 4; i++){
    for(int k = 0; k < 4; k++){
      if((bases[i] == 'T' && bases[k] == 'A') || (bases[i] == 'A' && bases[k] == 'T')){
        continue;
      }
      if((bases[i] == 'G' && bases[k] == 'C') || (bases[i] == 'C' && bases[k] == 'G')){
        continue;
      }
      print_num('G', 3);
      cout << char(bases[i]);
      print_num('G', 3);
      cout << "/";
      print_num('C', 3);
      cout << char(bases[k]);
      print_num('C', 3);
      cout << endl;
    }
  }

  for(int i = 0; i < 4; i++){
    for(int k = 0; k < 4; k++){
      if((bases[i] == 'T' && bases[k] == 'A') || (bases[i] == 'A' && bases[k] == 'T')){
        continue;
      }
      if((bases[i] == 'G' && bases[k] == 'C') || (bases[i] == 'C' && bases[k] == 'G')){
        continue;
      }
      print_num('A', 3);
      cout << char(bases[i]);
      print_num('A', 3);
      cout << "/";
      print_num('T', 3);
      cout << char(bases[k]);
      print_num('T', 3);
      cout << endl;
    }
  }
}
