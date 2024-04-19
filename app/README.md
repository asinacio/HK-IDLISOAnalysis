Directory to store source files for TreeConverter and CMake version of opticalAnalysis. A temporary version of analysis.cc is currently kept in the upper level for testing purposes, but on commit should be copied here and modified for compilation. Right now this should just involve replacing:
```
void runFit( std::string inFileName ){
```
with:
```
int main(int argc, char* argv[]){
```
