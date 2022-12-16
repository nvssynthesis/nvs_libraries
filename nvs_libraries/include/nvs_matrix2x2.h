//
//  nvs_matrix2x2.h
//  buttworth-efficiency-test
//
//  Created by Nicholas Solem on 12/13/22.
//  Copyright © 2022 Nicholas Solem. All rights reserved.
//

#ifndef nvs_matrix2x2_h
#define nvs_matrix2x2_h
#include <iostream>

namespace nvs_matrix {

struct vec2{
    float a, b;
    float *vec[2];

    vec2(){
        a = b = 0.f;
        vec[0] = &a;
        vec[1] = &b;
    }
    void print(){
        std::cout << a << "\n" << b << "\n";
    }

    static vec2 add(const vec2& v,const vec2& w){
        vec2 u;
        u.a = v.a + w.a;
        u.b = v.b + w.b;
        
        return u;
    }
    static vec2 scale(const vec2& v, const float&b){
        vec2 u;
        u.a = v.a * b;
        u.b = v.b * b;
        
        return u;
    }
    static float crossProduct(const vec2& v1, const vec2& v2){
        return (v1.a * v2.a) + (v1.b * v2.b);
    }
};

struct mat2x2{
    float a,b,c,d;
    float *mat[2][2];

    mat2x2(){
        a = b = c = d = 0.f;
        mat[0][0] = &a;
        mat[0][1] = &b;
        mat[1][0] = &c;
        mat[1][1] = &d;
    }
    
    void copy(const mat2x2& A){
        this->a = A.a;
        this->b = A.b;
        this->c = A.c;
        this->d = A.d;
    }
    void print(){
        for (int i = 0; i < 2; i++){
            for (int j = 0; j < 2; j++){
                std::cout << *(this->mat[i][j]) << "  ";
            }
            std::cout << "\n";
        }
    }
    
    static mat2x2 invert2x2(const mat2x2& A){
        mat2x2 B;
        float a = A.a;
        float b = A.b;
        float c = A.c;
        float d = A.d;
        
        float denom = a*d - b*c;
        if (denom == 0.f){
            std::cout << "NON INVERTIBLE\n";
            return B;
        }
        float scalar = 1.f / denom;
        
        B.a = scalar *  d;
        B.b = scalar * -b;
        B.c = scalar * -c;
        B.d = scalar *  a;

        return B;
    }
    
    static mat2x2 multiply2x2(const mat2x2& A, const mat2x2& B){
        mat2x2 C;
        
        C.a = A.a*B.a + A.b*B.c;
        C.b = A.a*B.b + A.b*B.d;
        C.c = A.c*B.a + A.d*B.c;
        C.d = A.c*B.b + A.d*B.d;
        
        return C;
    }
    static mat2x2 add2x2(const mat2x2& A, const mat2x2& B){
        mat2x2 C;
        
        C.a = A.a + B.a;
        C.b = A.b + B.b;
        C.c = A.c + B.c;
        C.d = A.d + B.d;
        
        return C;
    }
    static mat2x2 scale2x2(const mat2x2& A, const float& b){
        mat2x2 C;
        
        C.a = A.a * b;
        C.b = A.b * b;
        C.c = A.c * b;
        C.d = A.d * b;
        
        return C;
    }
    static vec2 matXvec(const mat2x2& X, const vec2& v){
        vec2 w;
        
        w.a = X.a*v.a + X.b*v.b;
        w.b = X.c*v.a + X.d*v.b;
        
        return w;
    }
};

} // nvs_matrix

#endif /* nvs_matrix2x2_h */
