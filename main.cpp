
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cctype>
#include <map>
#include <set>

using namespace std;

class Term {
public:
    int a, b, c, d;  // coefficient, x^b, sin^c x, cos^d x
    
    Term(int a = 0, int b = 0, int c = 0, int d = 0) : a(a), b(b), c(c), d(d) {}
    
    bool operator==(const Term& other) const {
        return b == other.b && c == other.c && d == other.d;
    }
    
    bool operator!=(const Term& other) const {
        return !(*this == other);
    }
    
    bool operator<(const Term& other) const {
        if (b != other.b) return b > other.b;
        if (c != other.c) return c > other.c;
        return d > other.d;
    }
    
    Term operator-() const {
        return Term(-a, b, c, d);
    }
};

class Poly {
public:
    vector<Term> terms;
    
    Poly() {}
    
    Poly(const vector<Term>& terms) : terms(terms) {
        simplify();
    }
    
    void simplify() {
        if (terms.empty()) return;
        
        // Sort terms by the required order
        sort(terms.begin(), terms.end());
        
        // Merge like terms
        vector<Term> result;
        for (const auto& term : terms) {
            if (term.a == 0) continue;
            
            if (result.empty() || result.back() != term) {
                result.push_back(term);
            } else {
                result.back().a += term.a;
            }
        }
        
        // Remove terms with zero coefficient
        terms.clear();
        for (const auto& term : result) {
            if (term.a != 0) {
                terms.push_back(term);
            }
        }
    }
    
    Poly operator+(const Poly& other) const {
        vector<Term> result = terms;
        result.insert(result.end(), other.terms.begin(), other.terms.end());
        return Poly(result);
    }
    
    Poly operator-(const Poly& other) const {
        vector<Term> negated;
        for (const auto& term : other.terms) {
            negated.push_back(-term);
        }
        vector<Term> result = terms;
        result.insert(result.end(), negated.begin(), negated.end());
        return Poly(result);
    }
    
    Poly operator*(const Poly& other) const {
        vector<Term> result;
        for (const auto& t1 : terms) {
            for (const auto& t2 : other.terms) {
                Term product;
                product.a = t1.a * t2.a;
                product.b = t1.b + t2.b;
                product.c = t1.c + t2.c;
                product.d = t1.d + t2.d;
                result.push_back(product);
            }
        }
        return Poly(result);
    }
    
    Poly derivate() const {
        vector<Term> result;
        for (const auto& term : terms) {
            if (term.a == 0) continue;
            
            // Derivative of x^b: b * x^(b-1)
            if (term.b > 0) {
                result.push_back(Term(term.a * term.b, term.b - 1, term.c, term.d));
            }
            
            // Derivative of sin^c x: c * sin^(c-1) x * cos x
            if (term.c > 0) {
                result.push_back(Term(term.a * term.c, term.b, term.c - 1, term.d + 1));
            }
            
            // Derivative of cos^d x: -d * sin x * cos^(d-1) x
            if (term.d > 0) {
                result.push_back(Term(-term.a * term.d, term.b, term.c + 1, term.d - 1));
            }
        }
        return Poly(result);
    }
    
    string toString() const {
        if (terms.empty()) return "0";
        
        string result;
        bool first = true;
        
        for (const auto& term : terms) {
            if (term.a == 0) continue;
            
            // Add operator
            if (!first) {
                if (term.a > 0) result += "+";
            }
            first = false;
            
            // Handle coefficient
            if (abs(term.a) == 1 && (term.b > 0 || term.c > 0 || term.d > 0)) {
                if (term.a < 0) result += "-";
            } else {
                result += to_string(term.a);
            }
            
            // Handle x^b
            if (term.b > 0) {
                result += "x";
                if (term.b > 1) {
                    result += "^" + to_string(term.b);
                }
            }
            
            // Handle sin^c x
            if (term.c > 0) {
                result += "sin";
                if (term.c > 1) {
                    result += "^" + to_string(term.c);
                }
                result += "x";
            }
            
            // Handle cos^d x
            if (term.d > 0) {
                result += "cos";
                if (term.d > 1) {
                    result += "^" + to_string(term.d);
                }
                result += "x";
            }
            
            // If the term is just a constant and it's negative, we need to handle the sign
            if (term.b == 0 && term.c == 0 && term.d == 0 && term.a < 0 && !result.empty() && result.back() == '-') {
                result.pop_back();
                result += to_string(term.a);
            }
        }
        
        return result;
    }
};

class Frac {
public:
    Poly p, q;  // p/q
    
    Frac() {}
    
    Frac(const Poly& p, const Poly& q) : p(p), q(q) {
        // Ensure denominator is not empty
        if (q.terms.empty()) {
            this->q = Poly({Term(1, 0, 0, 0)});
        }
    }
    
    Frac(int a) {
        p = Poly({Term(a, 0, 0, 0)});
        q = Poly({Term(1, 0, 0, 0)});
    }
    
    Frac(const Term& term) {
        p = Poly({term});
        q = Poly({Term(1, 0, 0, 0)});
    }
    
    Frac operator+(const Frac& other) const {
        Poly new_p = p * other.q + other.p * q;
        Poly new_q = q * other.q;
        return Frac(new_p, new_q);
    }
    
    Frac operator-(const Frac& other) const {
        Poly new_p = p * other.q - other.p * q;
        Poly new_q = q * other.q;
        return Frac(new_p, new_q);
    }
    
    Frac operator*(const Frac& other) const {
        Poly new_p = p * other.p;
        Poly new_q = q * other.q;
        return Frac(new_p, new_q);
    }
    
    Frac operator/(const Frac& other) const {
        Poly new_p = p * other.q;
        Poly new_q = q * other.p;
        return Frac(new_p, new_q);
    }
    
    Frac derivate() const {
        Poly p_prime = p.derivate();
        Poly q_prime = q.derivate();
        Poly num = p_prime * q - q_prime * p;
        Poly den = q * q;
        return Frac(num, den);
    }
    
    string toString() const {
        string p_str = p.toString();
        string q_str = q.toString();
        
        // Special cases
        if (p_str == "0") return "0";
        if (q_str == "1") return p_str;
        
        // Add parentheses if needed
        bool p_need_paren = p.terms.size() > 1;
        bool q_need_paren = q.terms.size() > 1;
        
        string result;
        if (p_need_paren) result += "(";
        result += p_str;
        if (p_need_paren) result += ")";
        
        result += "/";
        
        if (q_need_paren) result += "(";
        result += q_str;
        if (q_need_paren) result += ")";
        
        return result;
    }
};

// Parser for the expression
class Parser {
private:
    string s;
    int pos;
    
    // Skip whitespace
    void skipSpace() {
        while (pos < s.length() && isspace(s[pos])) {
            pos++;
        }
    }
    
    // Get number from current position
    int getNumber() {
        skipSpace();
        if (pos >= s.length()) return 1;
        
        int sign = 1;
        if (s[pos] == '-') {
            sign = -1;
            pos++;
            skipSpace();
        } else if (s[pos] == '+') {
            pos++;
            skipSpace();
        }
        
        // Check if there's a digit
        if (pos < s.length() && isdigit(s[pos])) {
            int num = 0;
            while (pos < s.length() && isdigit(s[pos])) {
                num = num * 10 + (s[pos] - '0');
                pos++;
            }
            return sign * num;
        }
        
        // No digit, return ±1
        return sign;
    }
    
    // Parse exponent
    int parseExp() {
        skipSpace();
        if (pos < s.length() && s[pos] == '^') {
            pos++;
            skipSpace();
            int num = 0;
            while (pos < s.length() && isdigit(s[pos])) {
                num = num * 10 + (s[pos] - '0');
                pos++;
            }
            return num;
        }
        return 1;
    }
    
    // Parse a term
    Term parseTerm() {
        skipSpace();
        if (pos >= s.length()) return Term(1, 0, 0, 0);
        
        int coeff = getNumber();
        
        int x_exp = 0, sin_exp = 0, cos_exp = 0;
        
        // Parse x
        skipSpace();
        if (pos < s.length() && s[pos] == 'x') {
            pos++;
            x_exp = parseExp();
            if (x_exp == 0) x_exp = 1;  // x without ^ means x^1
        }
        
        // Parse sin
        skipSpace();
        if (pos + 2 < s.length() && s.substr(pos, 3) == "sin") {
            pos += 3;
            sin_exp = parseExp();
            if (sin_exp == 0) sin_exp = 1;  // sin without ^ means sin^1
            skipSpace();
            if (pos < s.length() && s[pos] == 'x') {
                pos++;
            }
        }
        
        // Parse cos
        skipSpace();
        if (pos + 2 < s.length() && s.substr(pos, 3) == "cos") {
            pos += 3;
            cos_exp = parseExp();
            if (cos_exp == 0) cos_exp = 1;  // cos without ^ means cos^1
            skipSpace();
            if (pos < s.length() && s[pos] == 'x') {
                pos++;
            }
        }
        
        return Term(coeff, x_exp, sin_exp, cos_exp);
    }
    
    // Parse a factor (term or parenthesized expression)
    Frac parseFactor() {
        skipSpace();
        if (pos >= s.length()) return Frac(1);
        
        if (s[pos] == '(') {
            pos++;  // skip '('
            Frac result = parseExpr();
            skipSpace();
            if (pos < s.length() && s[pos] == ')') {
                pos++;  // skip ')'
            }
            return result;
        }
        
        return Frac(parseTerm());
    }
    
    // Parse multiplication and division
    Frac parseMulDiv() {
        Frac result = parseFactor();
        
        while (true) {
            skipSpace();
            if (pos < s.length() && s[pos] == '*') {
                pos++;
                Frac rhs = parseFactor();
                result = result * rhs;
            } else if (pos < s.length() && s[pos] == '/') {
                pos++;
                Frac rhs = parseFactor();
                result = result / rhs;
            } else {
                break;
            }
        }
        
        return result;
    }
    
    // Parse addition and subtraction
    Frac parseExpr() {
        Frac result = parseMulDiv();
        
        while (true) {
            skipSpace();
            if (pos < s.length() && s[pos] == '+') {
                pos++;
                Frac rhs = parseMulDiv();
                result = result + rhs;
            } else if (pos < s.length() && s[pos] == '-') {
                pos++;
                Frac rhs = parseMulDiv();
                result = result - rhs;
            } else {
                break;
            }
        }
        
        return result;
    }
    
public:
    Frac parse(const string& expr) {
        s = expr;
        pos = 0;
        return parseExpr();
    }
};

void solve(const string& expr) {
    Parser parser;
    Frac f = parser.parse(expr);
    
    // Output the original expression as a fraction
    cout << f.toString() << endl;
    
    // Output the derivative
    Frac f_prime = f.derivate();
    cout << f_prime.toString() << endl;
}

int main() {
    string expr;
    getline(cin, expr);
    solve(expr);
    return 0;
}
