#include <iostream>
#include <vector>
#include <cmath>
using namespace std;
vector<vector<double>> E(int n){ //identity matrix of size n
    vector<vector<double>> result (n);
    for (int i = 0; i<n; i++){
        result[i].resize(n,0);
        result[i][i] = 1;
    }
    return result;
}
vector <vector<double>> orth(int n){ //unit vector of size n
    vector <vector<double>> result(n);
    result[0].push_back(1);
    for(int i = 1; i < n; i++){
        result[i].push_back(0);
    }
    return result;
}
vector<double> getcolumn (vector<vector<double>> A, int k){ //get a column of number k as a vector from a matrix
    vector<double> result;
    for(int i = 0; i < A.size(); i++){
        result.push_back(A[i][k]);
    }
    return result;
}
vector <vector<double>> msub (vector <vector<double>> A, vector <vector<double>> B){ //subtract matrix B from matrix A
    vector<vector<double>> result (A.size());
    for (int i = 0; i < result.size(); i++){
        result[i].resize(A[0].size());
        for (int j = 0; j < result[i].size(); j++){
            result[i][j] = A[i][j]-B[i][j];
        }
    }
    return result;
}
double ger_sqrt(double u, double eps){ //square root
    double res = 1, old = 0;
    int k = 0;
    while (abs(res - old) > eps){
        old = res;
        res = 0.5*(old+double(u)/old);
    }
    return res;
}
vector <vector<double>> mnumm(double a, vector <vector<double>> b){ //multiply a matrix by a number
    for (int i = 0; i < b.size(); i++){
        for(int j = 0; j < b[i].size(); j++) {
            b[i][j] = b[i][j]*a;
        }
    }
    return b;
}
vector <vector<double>> mmultiply (vector <vector<double>> A, vector <vector<double>> B){ //matrix multiplication
    if (A[0].size() == B.size()){
        vector <vector<double>> C (A.size());
        for (int i = 0; i < C.size(); i++){
            C[i].resize(B[0].size());
        }
        for (int i = 0; i < A.size(); i++){
            for (int j = 0; j < B[0].size(); j++){
                for (int k = 0; k < A[0].size(); k++){
                    C[i][j] += A[i][k]*B[k][j];
                }
            }
        }
        return C;
    }
}
double vnorm (vector<vector<double>> a){ //euclidian norm of a vector
    double result = 0;
    for (auto i:a){
        result+=i[0]*i[0];
    }
    return ger_sqrt(result, pow(10,-10));
}
void normalize (vector<vector<double>>& a){ //get a unit vector of the vector
    double p = vnorm(a);
    for (int i = 0; i < a.size(); i++){
        a[i][0] /= p;
    }
}
double NormInf(vector<vector<double>> A) { //  maximum absolute row sum of the matrix.
    double result = INT64_MIN, s = 0;
    for (int i = 0; i < A.size(); i++){
        for (int j = 0; j < A.size(); j++){
            s += abs(A[i][j]);
        }
        if (s > result){
            result = s;
            s = INT64_MIN;
        }
        else{
            s = INT64_MIN;
        }
    }
    return result;
}
double NormInfCol(vector<vector<double>> A) { //maximum absolute column sum of the matrix
    double result = INT64_MIN, s = 0;
    for (int i = 0; i < A.size(); i++){
        for (int j = 0; j < A.size(); j++){
            s += abs(A[j][i]);
        }
        if (s > result){
            result = s;
            s = INT64_MIN;
        }
        else{
            s = INT64_MIN;
        }
    }
    return result;
}
vector<vector<double>> transpose (vector<vector<double>> A){ //returns the transposed matrix
    vector<vector<double>> result (A[0].size());
    for(int i = 0; i < result.size(); i++){
        result[i] = getcolumn(A,i);
    }
    return result;
}
vector<vector<double>> getblock(vector<vector<double>> A, int a, int b){ //get a square submatrix of the matrix from the element of the (a,a) position to the element of (b,b) position
    vector<vector<double>> result(b-a+1);
    for(int i = 0; i < b-a+1; i++){
        result[i].resize(b-a+1);
        for (int j = 0; j < b-a+1; j++){
            result[i][j] = A[i+a][j+a];
        }
    }
    return result;
}
void insertblock(vector<vector<double>>& A, vector<vector<double>>B, int a){ //insert a block into your matrix with the (1,1) element of the block going into (a,a) element of matrix A
    for (int i = a; i < a + B.size(); i++){
        for (int j = a; j < a + B.size(); j++){
            A[i][j] = B[i-a][j-a];
        }
    }
}
pair<vector<vector<double>>,vector<vector<double>>> LUP (vector <vector<double>> A){ //returns a matrix that stores L,U and a matrix P (transposition matrix)
    vector <vector<double>> P(A.size());
    for(int i = 0; i < P.size(); i++){
        P[i].resize(A.size(),0);
        P[i][i] = 1;
    }
    for(int i = 0; i<A.size()-1; i++){
        double lead = INT64_MIN;
        double nlead = -1;
        for (int j = i; j < A.size(); j++){
            if (abs(A[j][i]) > lead){
                lead = abs(A[j][i]);
                nlead = j;
            }
        }
        swap(A[i],A[nlead]);
        swap(P[i],P[nlead]);
        for (int j = i+1; j < A.size(); j++){
            A[j][i] = A[j][i]/A[i][i];
            for (int k = i+1; k<A.size(); k++){
                A[j][k] = A[j][k]-A[j][i]*A[i][k];
            }
        }
    }
    return make_pair(A,P);
}
vector <double> LUPsolve(vector <vector<double>> A, vector<vector<double>> b){ //solves the equation system by using the results of LUP function
    pair<vector<vector<double>>,vector<vector<double>>> LpUaP = LUP (A);
    vector<vector<double>> LU = LpUaP.first;
    b = mmultiply(LpUaP.second,b);
    vector<double> y(b.size());
    for(int i = 0; i<b.size(); i++){
        y[i] = b[i][0];
    }
    for (int i = 0; i < A.size(); i++){
        for (int k = 0; k<i;k++){
            y[i]-=LU[i][k]*y[k];
        }
    }
    vector<double> x(b.size());
    for(int i = b.size()-1; i>=0; i--){
        x[i] = y[i];
        for (int k = i+1; k<b.size(); k++){
            x[i] -= LU[i][k]*x[k];
        }
        x[i] = x[i]/LU[i][i];
    }
    return x;
}
pair <vector <vector<double>>, vector <vector<double>>> QR (vector <vector<double>> A, vector <vector<double>> b){ //Returns Q and R matrixes
    int n = A.size();
    vector <vector<double>> Q = E(n), R = A, w, z, y;
    for (int i = 0; i < n-1; i++){
        y.resize(n-i);
        for(int j = 0; j < y.size(); j++){
            y[j].push_back(R[i+j][i]);
        }
        double a = vnorm(y);
        z = orth(n-i);
        vector<vector<double>> Qi = E(n), Ri;
        w = msub(y, mnumm(a,z));
        normalize(w);
        insertblock(Qi, msub(E(Q.size()-i),mnumm(2, mmultiply(w, transpose(w)))),i);
        Q = mmultiply(Q,Qi);
        Ri = mmultiply(getblock(Qi,i,2), getblock(R,i,2));
        insertblock(R,Ri,i);
        y.resize(0);
    }
    return make_pair(Q,R);
}
vector<double> QRsolve (vector <vector<double>> A, vector <vector<double>> b){ //returns a solution using the results of QR function
    pair <vector <vector<double>>, vector <vector<double>>> decomp = QR(A,b);
    vector<double> y = getcolumn(mmultiply(transpose(decomp.first),b),0);
    vector<double> x = y;
    for(int i = b.size()-1; i>=0; i--){
        for (int k = i+1; k<b.size(); k++){
            x[i] -= decomp.second[i][k]*x[k];
        }
        x[i] = x[i]/decomp.second[i][i];
    }
    return x;
}
vector<double> FPISolve(vector <vector<double>> A, vector<vector<double>> b, double eps){ //returns a solution approximated using fixed-point iterations
    double mu = 1/ NormInfCol(A);
    vector<vector<double>> B = msub(E(A.size()),mnumm(mu,A));
    double (*mnorm)(vector<vector<double>>) = nullptr;
    if(NormInf(B) < 1){
        mnorm = &NormInf;
    }
    else{
        if (NormInfCol(B) < 1){
            mnorm = &NormInfCol;
        }
        else{
            vector<vector<double>> AT = transpose(A);
            A = mmultiply(AT,A);
            b = mmultiply(AT,b);
            AT.clear();
            mu = 1/ NormInfCol(A);
            B = msub(E(A.size()),mnumm(mu,A));
            if(NormInf(B) < 1){
                mnorm = &NormInf;
            }
            else {
                if (NormInfCol(B) < 1) {
                    mnorm = &NormInfCol;
                }
            }
        }
    }
    vector<vector<double>> c = mnumm(mu,b);
    vector<vector<double>> x, xprev;
    x = c;
    double k = 0;
    if(mnorm!= nullptr){
        double normb = mnorm(B);
        do{
            xprev = x;
            x = msub(mmultiply(B,x), mnumm(-1,c));
            k++;
        }while((normb/(1-normb))* vnorm(msub(x,xprev)) > eps);
    }
    else{
        do{
            x = msub(mmultiply(B,x), mnumm(-1,c));
            k++;
        }while(vnorm(msub(mmultiply(A,x),b)) > eps);
    }
    x.push_back({k});
    return getcolumn(x,0);
}
bool ddom (vector<vector<double>> A){ //checks if the matrix is diagonally dominant
    for(int i = 0; i < A.size(); i++){
        double c = 0;
        for (int j = 0; j < A.size(); j++){
            c+=abs(A[i][j]);
        }
        if(2*abs(A[i][i]) < c){
            return false;
        }
    }
    return true;
}
vector <double> SeidelSolve (vector <vector<double>> A, vector<vector<double>> b, double eps){ //returns a solution approximated using Gauss-Seidel iterative process.
    if (!ddom(A)){
        vector<vector<double>> AT = transpose(A);
        A = mmultiply(AT,A);
        b = mmultiply(AT,b);
        AT.clear();
    }
    vector<vector<double>>x;
    vector <vector<double>> C = A;
    vector <vector<double>> d = b;
    for(int i = 0; i < C.size(); i++){
        double diag = C[i][i];
        C[i][i] = 0;
        for(int j = 0; j <C.size(); j++){
            if(i!=j){
                C[i][j] = -C[i][j]/diag;
            }
        }
        d[i][0] = d[i][0]/diag;
    }
    x = d;
    double k = 0;
    do {
        for (int i = 0; i < C.size(); i++) {
            double xtemp = 0;
            for (int j = 0; j < C.size(); j++) {
                xtemp += C[i][j] * x[j][0];
            }
            x[i][0] = xtemp;
            x[i][0] += d[i][0];
        }
        k++;
    }while(vnorm(msub(mmultiply(A,x),b)) > eps);
    x.push_back({k});
    return getcolumn(x,0);
}
int main() {
    int n;
    cin >> n;
    vector <vector<double>> A(n);
    vector <vector<double>> b(n);
    for (int i = 0; i < n; i++){
        A[i].resize(n);
        b[i].resize(1);
        cin >> b[i][0];
    }
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            cin >> A[i][j];
        }
    }
    cout << "LU decomposition" << endl;
    vector <double> x = LUPsolve(A,b);
    for (int i = 0; i < n; i++){
        cout << "x" << i << " = " << x[i] << endl;
    }
    cout << "QR decomposition" << endl;
    x = QRsolve(A,b);
    for (int i = 0; i < n; i++){
        cout << "x" << i << " = " << x[i] << endl;
    }
    x = FPISolve(A,b,pow(10,-3));
    cout << "Fixed-point iterations" << endl;
    for (int i = 0; i < n; i++){
        cout << "x" << i << " = " << x[i] << endl;
    }
    cout << x[n] << " iterations" << endl;
    x = SeidelSolve(A,b,pow(10,-3));
    for (int i = 0; i < n; i++){
        cout << "x" << i << " = " << x[i] << endl;
    }
    cout << x[n] << " iterations" << endl;
}
