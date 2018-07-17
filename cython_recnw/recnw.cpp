
/*-------------------------------------------
 * AUTHOR: Alexandre YAHI
 * AFFILIATION: Columbia University Medical Center
 * 2015-2018
 -------------------------------------------*/

#include "recnw.h"
#include <utility>
#include <map>
#include <sstream>
#include <limits>

using namespace std;

const float minus_inf = - std::numeric_limits<float>::infinity();

void  F_init_aff( float ** F, int L1, int L2, float gap_op, float gap_ext, bool free_hgap_1, bool free_hgap_2)
{ // Affine gap penalty initalization for the match matrix
        if(free_hgap_1 || free_hgap_2){ // at least one free head gap
            F[ 0 ][ 0 ] =  0.0 ;
        }else{
            F[ 0 ][ 0 ] =  minus_inf; // no free head gap
        }

        int i=0, j=0;

        if(free_hgap_1){
            for( j = 1; j < L2; j++ )
            {
                    F[ 0 ][ j ] =  0.0 ; // Free head gap in seq1
            }
        }else{
            for( j = 1; j < L2; j++ )
            {
                    // F[ 0 ][ j ] =  -32767 ; // Standard
                    F[ 0 ][ j ] =  - (gap_op + j*gap_ext) ; // Standard
            }
        }

        if(free_hgap_2){
            for( i = 1; i < L1; i++ )
            {
                    F[ i ][ 0 ] =  0.0 ; // Free head gap in seq2
            }
        }else{
            for( i = 1; i < L1; i++ )
            {
                    // F[ i ][ 0 ] =  -32767 ; // Standard
                    F[ i ][ 0 ] =  - (gap_op + i*gap_ext) ; // Standard
            }
        }
}


void  F_init_lin( float ** F, int L1, int L2, float gap_penalty, bool free_hgap_1, bool free_hgap_2)
{ // Linear gap penalty - free head gaps
        F[ 0 ][ 0 ] =  0.0 ;

        int i=0, j=0;

        if(free_hgap_1){
            for( j = 1; j < L2; j++ )
            {
                    F[ 0 ][ j ] = -j*gap_penalty ;
            }
        }else{
            for( j = 1; j < L2; j++ )
            {
                    F[ 0 ][ j ] = 0.0 ;
            }
        }

        if(free_hgap_2){
            for( i = 1; i < L1; i++ )
            {
                    F[ i ][ 0 ] = 0.0 ;
            }
        }else{
            for( i = 1; i < L1; i++ )
            {
                    F[ i ][ 0 ] = -i*gap_penalty ;
            }
        }
}


void  P_init( float ** P, int L1, int L2, float gap_op, float gap_ext, bool free_hgap_2)
{ // Gap in seq2 matrix initialization
        P[ 0 ][ 0 ] =  minus_inf ;

        int i=0, j=0;

        // First row
        for( j = 1; j < L2; j++ )
        {
                P[ 0 ][ j ] =  minus_inf ;
        }

        // First column
        if(free_hgap_2){
            for( i = 1; i < L1; i++ )
            {
                    P[ i ][ 0 ] =  0.0; // Free head gap in seq2
            }
        }else{
            for( i = 1; i < L1; i++ )
            {
                    P[ i ][ 0 ] =  - (gap_op + i*gap_ext) ; // Standard
            }
        }

}

void  Q_init( float ** Q, int L1, int L2, float gap_op, float gap_ext, bool free_hgap_1)
{ // Gap in seq1 matrix initialization
        Q[ 0 ][ 0 ] =  minus_inf ;

        int i=0, j=0;

        // First row
        if(free_hgap_1){
            for( j = 1; j < L2; j++ )
            {
                    Q[ 0 ][ j ] =  0.0 ; // Free head gap in seq1
            }
        }else{
            for( j = 1; j < L2; j++ )
            {
                    Q[ 0 ][ j ] =  - (gap_op + j*gap_ext) ; // Standard
            }
        }

        // First column
        for( i = 1; i < L1; i++ )
        {
                Q[ i ][ 0 ] =  minus_inf ; // Semi-global variant
        }
}


float max_2(float a, float b)
{
    float max = 0.0;
    if(a >= b){
        max = a;
    }
    else{
        max = b;
    }
    return max;
}


float max_3(float fU, float fD, float fL)
{
    float max = 0.0;
    max = max_2(fD, fL) ;
    max = max_2(fU, max) ;
    return max;
}

/// MASTER FUNCTION - AFFINE
std::string recnw_affine(
              string seq_1,            /* fixed sequence - seq1 - reference */
              string seq_2,              /* read to align - seq2 - read */
              float gap_op,                   /* gap opening penality, positive int */
              float gap_ext,                   /* gap extension penalty, positive int */
              float match,                    /* Match score, positive int */
              float mismatch,                  /* Mismatch penalty, negative int*/
              bool free_hgap_1,               /* No head gap penalty in seq1 */
              bool free_hgap_2,               /* No head gap penalty in seq2 */
              bool free_tgap_1,               /* No tail gap penalty in seq1 */
              bool free_tgap_2,               /* No tail gap penalty in seq2 */
              int sim,                      /* Similarity with previous run */
              int terminate                 /* clear memory */
            )
    {
        // // CONST AND VAR
        const int  L1 = seq_1.length()+1;
        const int  L2 = seq_2.length()+1;
        float        dg_op, dg_ext, ig_op, ig_ext, no_gap ;
        int        i = 0, j = 0;

        // Dynamic programming matrix - INITALIZATION
        static float ** F;
        static float ** P;
        static float ** Q;

        if(sim==-1){             // Only initialize if first round
            F = new float * [L1];
            P = new float * [L1];
            Q = new float * [L1];
            for( int i = 0; i < L1; i++ ){
                F[ i ] = new float[L2];
                P[ i ] = new float[L2];
                Q[ i ] = new float[L2];
            }

            F_init_aff(F, L1, L2, gap_op, gap_ext, free_hgap_1, free_hgap_2);
            P_init(P, L1, L2, gap_op, gap_ext, free_hgap_2);
            Q_init(Q, L1, L2, gap_op, gap_ext, free_hgap_1);
        }

        i=0, j=0;

        string seq_1_al, seq_2_al;

        std::map <char,int> m; m['A']=0; m['C']=1; m['G']=2;m['T']=3;m['N']=4;

        const float  a =  match;   // Match
        const float  b = mismatch;   // Mismatch

        const float  s[ 5 ][ 5 ] = { { a, b, b, b, -2.0 },    /* substitution matrix */
                                   { b, a, b, b, -2.0 },
                                   { b, b, a, b, -2.0 },
                                   { b, b, b, a, -2.0 },
                                   {-2.0,-2.0,-2.0,-2.0, -1.0 }} ;

        if(sim==-1){
            for( i = 1; i < L1; i++ )
            {
                for( j = 1; j < L2; j++ )
                {
                    // Take maximum because penalties are negative
                    //Compute first P the deletion gap-opening/gap-extension matrix
                    dg_op = F[i-1][j] - (gap_op);
                    // dg_op = F[i-1][j] - (gap_op + gap_ext);
                    dg_ext = P[i-1][j] - gap_ext;

                    P[i][j] = max_2(dg_op, dg_ext); // DELETION

                    ig_op = F[i][j-1] - (gap_op);
                    // ig_op = F[i][j-1] - (gap_op + gap_ext);
                    ig_ext = Q[i][j-1] - gap_ext;

                    Q[i][j] = max_2(ig_op, ig_ext); // INSERTION

                    no_gap = F[ i-1 ][ j-1 ] + s[m[seq_1[i-1]]][m[seq_2[j-1]]]; // MATCH

                    F[ i ][ j ] = max_3(no_gap, P[i][j], Q[i][j]) ;
                }
            }
        }
        else{
           for( i = 1; i < L1; i++ )
           {
               for( j = sim + 1; j < L2; j++ )
               {
                   // Take maximum because penalties are negative
                   //Compute first P the deletion gap-opening/gap-extension matrix
                   dg_op = F[i-1][j] - (gap_op + gap_ext);
                   dg_ext = P[i-1][j] - gap_ext;

                   P[i][j] = max_2(dg_op, dg_ext); // DELETION

                   ig_op = F[i][j-1] - (gap_op + gap_ext);
                   ig_ext = Q[i][j-1] - gap_ext;

                   Q[i][j] = max_2(ig_op, ig_ext); // INSERTION

                   no_gap = F[ i-1 ][ j-1 ] + s[m[seq_1[i-1]]][m[seq_2[j-1]]]; // MATCH

                   F[ i ][ j ] = max_3(no_gap, P[i][j], Q[i][j]) ;
               }
           }
       }

       i-- ; j-- ;

        // TRACEBACK

        int i_max = 0; // index i of best score
        int fi_max = 0; // best score finishing by a match
        int pi_max = 0; // best score finishing in a gap in seq2

        int j_max = 0;
        int qj_max = 0;
        int fj_max = 0;

        float temp_i = 0;
        float temp_fi = 0;
        float temp_pi = 0;

        float temp_j = 0;
        float temp_fj = 0;
        float temp_qj = 0;

        float score = 0;
        float score_i = 0;
        float score_j = 0;

        i = L1-1;
        j = L2-1;


        if(free_tgap_1 && free_tgap_2 == false){ // Free trailing in seq 1
            // Look for max in last row - Free tail in seq1
            temp_qj = minus_inf;
            for(int l = 1; l < L2; l++)
            {
                if( Q[L1-1][l] > temp_qj )
                {
                    qj_max = l;
                    temp_qj = Q[L1-1][l];
                }
            }

            temp_fj = minus_inf;
            for(int l = 1; l < L2; l++)
            {
                if( F[L1-1][l] > temp_fj )
                {
                    fj_max = l;
                    temp_fj = F[L1-1][l];
                }
            }

            // Get max between Q and F last row - free trail seq1
            if(temp_qj > temp_fj){
                j_max = qj_max;
                score = Q[L1-1][j_max];
            }else{
                j_max = fj_max;
                score = F[L1-1][j_max];
            }

            // Initiate the traceback
            i_max = L1-1;
            while(j>j_max)
            {
              seq_1_al = '-' + seq_1_al;
              seq_2_al = seq_2[j-1] + seq_2_al;
              j--;
            }

        }else if(free_tgap_2 && free_tgap_1 == false){ // Free trailing in seq 2
            // Look for max in the last column - Free tail in seq2
            temp_pi = minus_inf;
            for(int k = 1; k < L1; k++)
            {
                if( P[k][L2-1] >= temp_pi )
                {
                    pi_max = k;
                    temp_pi = P[k][L2-1];
                }
            }

            temp_fi = minus_inf;
            for(int k = 1; k < L1; k++)
            {
                if( F[k][L2-1] >= temp_fi )
                {
                    fi_max = k;
                    temp_fi = F[k][L2-1];
                }
            }

            // Get max between P and F last column - free trail seq2
            j_max = L2-1;
            if(temp_pi > temp_fi){
                i_max = pi_max;
                score = P[i_max][L2-1];
            }else{
                i_max = fi_max;
                score = F[i_max][L2-1];
            }

            // Initiate the traceback
            j_max = L2-1;
            while(i>i_max)
            {
                seq_1_al = seq_1[i-1] + seq_1_al;
                seq_2_al = '-' + seq_2_al;
                i--;
            }
        }else if(free_tgap_1 && free_tgap_2){ // Free trailing in both
            // Look for max in the last column - Free tail in seq2
            temp_pi = minus_inf;
            for(int k = 1; k < L1; k++)
            {
                if( P[k][L2-1] >= temp_pi )
                {
                    pi_max = k;
                    temp_pi = P[k][L2-1];
                }
            }

            temp_fi = minus_inf;
            for(int k = 1; k < L1; k++)
            {
                if( F[k][L2-1] >= temp_fi )
                {
                    fi_max = k;
                    temp_fi = F[k][L2-1];
                }
            }

            // Look for max in last row - Free tail in seq1
            temp_qj = minus_inf;
            for(int l = 1; l < L2; l++)
            {
                if( Q[L1-1][l] > temp_qj )
                {
                    qj_max = l;
                    temp_qj = Q[L1-1][l];
                }
            }

            temp_fj = minus_inf;
            for(int l = 1; l < L2; l++)
            {
                if( F[L1-1][l] > temp_fj )
                {
                    fj_max = l;
                    temp_fj = F[L1-1][l];
                }
            }

            // Get max between P and F last column - free trail seq2
            j_max = L2-1;
            if(temp_pi > temp_fi){
                i_max = pi_max;
                score_i = P[i_max][L2-1];
                temp_i = temp_pi;
            }else{
                i_max = fi_max;
                score_i = F[i_max][L2-1];
                temp_i = temp_fi;
            }

            // Get max between Q and F last row - free trail seq1
            if(temp_qj > temp_fj){
                j_max = qj_max;
                score_j = Q[L1-1][j_max];
                temp_j = temp_qj;
            }else{
                j_max = fj_max;
                score_j = F[L1-1][j_max];
                temp_j = temp_fj;
            }

            // Decide between free trail in seq1 or seq2
            if(temp_i > temp_j){
                j_max = L2-1;
                score = score_i;
                while(i>i_max)
                {
                    seq_1_al = seq_1[i-1] + seq_1_al;
                    seq_2_al = '-' + seq_2_al;
                    i--;
                }
            }else{
                i_max = L1-1;
                score = score_j;
                while(j>j_max)
                {
                  seq_1_al = '-' + seq_1_al;
                  seq_2_al = seq_2[j-1] + seq_2_al;
                  j--;
                }
            }
        }else{ // No free trailing
            i_max = L1-1;
            j_max = L2-1;
            score = F[i_max][j_max];
        }


        int mat = 0; //0 is F, 1 is P, 2 is Q - help to jump between matrices

        // TRANSVERSAL TRACEBACK
        while( (i > 0) && (j > 0) )
        {
            if(mat == 0)
            {
                if(F[i][j] == F[i-1][j-1] + s[m[seq_1[i-1]]][m[seq_2[j-1]]])
                {
                    seq_1_al = seq_1[i-1] + seq_1_al;
                    seq_2_al = seq_2[j-1] + seq_2_al;
                    i--;
                    j--;
                }
                else if(F[i][j] == P[i][j])
                {
                    mat = 1; // Jump to matrix P for deletion
                }
                else
                {
                    mat = 2; // Jump to matrix Q for insertion
                }
            }
            else if(mat == 1) // DELETION - matrix P
            {
                if(P[i][j]==P[i-1][j] - gap_ext) //Extension, stay in P
                {
                    seq_1_al = seq_1[i-1] + seq_1_al;
                    seq_2_al = '-' + seq_2_al;
                    i--;
                    mat = 1;
                }
                else // Opening - move back to F
                {
                    seq_1_al = seq_1[i-1] + seq_1_al;
                    seq_2_al = '-' + seq_2_al;
                    i--;
                    mat = 0;
                }
            }
            else // INSERTION - matrix Q
            {
                if(Q[i][j] == Q[i][j-1] - gap_ext) //Extension - stay in Q
                {
                    seq_1_al = '-' + seq_1_al;
                    seq_2_al = seq_2[j-1] + seq_2_al;
                    j--;
                    mat = 2;
                }
                else //Opening - move back to F
                {
                    seq_1_al = '-' + seq_1_al;
                    seq_2_al = seq_2[j-1] + seq_2_al;
                    j--;
                    mat = 0;
                }
            }
        }


        if(j==0)
        {
            while(i>0)
            {
                seq_1_al = seq_1[i-1] + seq_1_al;
                seq_2_al = '-' + seq_2_al;
                i--;
            }
        }
        else
        {
            while(j>0)
            {
                seq_1_al = '-' + seq_1_al;
                seq_2_al = seq_2[j-1] + seq_2_al;
                j--;
            }
        }

        // formating of the score
        string Score;
        ostringstream convert;
        convert << score;
        Score = convert.str();


        string alignment = seq_1_al + '|' + seq_2_al + '|' + Score ;
        //De-allocate memory if terminate
        if (terminate == 1)
        {
            for (int i = 0; i < L1; ++i)
                {
                    delete [] F[i];
                    delete [] P[i];
                    delete [] Q[i];
                }

            delete [] F;
            delete [] P;
            delete [] Q;
        }

        return alignment;
}


/// MASTER FUNCTION - REGULAR
std::string recnw_reg(
              string seq_1,             /* fixed sequence */
              string seq_2,             /* read to align */
              float gap_penalty,                   /* gap penality, positive int */
              float match,                    /* Match score, positive int */
              float mismatch,                  /* Mismatch penalty, negative int*/
              bool free_hgap_1,               /* No head gap penalty in seq1 */
              bool free_hgap_2,               /* No head gap penalty in seq2 */
              bool free_tgap_1,               /* No tail gap penalty in seq1 */
              bool free_tgap_2,               /* No tail gap penalty in seq2 */
              int sim,                      /* Similarity with previous run */
              int terminate                 /* clear memory */
            )
    {
        // // CONST AND VAR
        const int  L1 = seq_1.length()+1;
        const int  L2 = seq_2.length()+1;
        float        del_cost, ins_cost, match_cost;
        int        i = 0, j = 0;

        // Dynamic programming matrix - INITALIZATION
        static float ** F;

        if(sim==-1){       // Only initialize if first round
            F = new float * [L1];
            for( int i = 0; i < L1; i++ ){
                F[ i ] = new float [L2];
            }

            F_init_lin(F, L1, L2, gap_penalty, free_hgap_1, free_hgap_2);
        }
        //

        i=0, j=0;

        string seq_1_al, seq_2_al;

        std::map <char,int> m; m['A']=0; m['C']=1; m['G']=2;m['T']=3;m['N']=4;

        const float  a = match;   // Match
        const float  b = mismatch;   // Mismatch

        const float  s[ 5 ][ 5 ] = { { a, b, b, b, -2.0 },    /* substitution matrix */
                                   { b, a, b, b, -2.0 },
                                   { b, b, a, b, -2.0 },
                                   { b, b, b, a, -2.0 },
                                   {-2.0,-2.0,-2.0,-2.0, -1.0 }} ;


        if(sim < 0){
            for( i = 1; i < L1; i++ )
            {
                for( j = 1; j < L2; j++ )
                {
                    del_cost = F[i-1][j] - gap_penalty ; // DELETION
                    match_cost = F[i-1][j-1] + s[m[seq_1[i-1]]][m[seq_2[j-1]]]; // MATCH
                    ins_cost = F[i][j-1] - gap_penalty ; // INSERTION

                    F[i][j] = max_3(del_cost, match_cost, ins_cost) ;
                }
            }
        }else{ // When re-using lines already computed
            for( i = 1; i < L1; i++ )
            {
                for( j = sim + 1; j < L2; j++ )
                {
                    del_cost = F[i-1][j] - gap_penalty ; // DELETION
                    match_cost = F[i-1][j-1] + s[m[seq_1[i-1]]][m[seq_2[j-1]]]; // MATCH
                    ins_cost = F[i][j-1] - gap_penalty ; // INSERTION

                    F[i][j] = max_3(del_cost, match_cost, ins_cost) ;
                }
            }
        }

        i-- ; j-- ;

        // FREE TRAILING TRACEBACK IMPLEMENTATION
        int i_max = 0;
        int j_max = 0;
        float temp_i = 0.0;
        float temp_j = 0.0;

        float score = 0.0;

        i = L1-1;
        j = L2-1;

        if(free_tgap_1 && free_tgap_2 == false){ // Free trailing seq 1
            // Look for max in the last row
            temp_j = minus_inf;
            for(int l = 1; l < L2; l++)
            {
                if( F[L1-1][l] > temp_j )
                {
                    j_max = l;
                    temp_j = F[L1-1][l];
                }
            }

            // Start adding gap in seq1 if needed
            i_max = L1-1;
            while(j>j_max)
            {
                seq_1_al = '-' + seq_1_al;
                seq_2_al = seq_2[j-1] + seq_2_al;
                j--;
            }
        }else if(free_tgap_2 && free_tgap_1 == false){ // Free trailing seq 2
            temp_i = minus_inf;
            for(int k = 1; k < L1; k++)
            {
                if( F[k][L2-1] > temp_i )
                {
                    i_max = k;
                    temp_i = F[k][L2-1];
                }
            }

            // Start adding gap in seq2 if needed
            j_max = L2-1;
            while(i>i_max)
            {
                seq_1_al = seq_1[i-1] + seq_1_al;
                seq_2_al = '-' + seq_2_al;
                i--;
            }

        }else if(free_tgap_1 && free_tgap_2){ // Free trailing in both
            temp_i = minus_inf;
            for(int k = 1; k < L1; k++)
            {
                if( F[k][L2-1] > temp_i )
                {
                    i_max = k;
                    temp_i = F[k][L2-1];
                }
            }

            // Look for max in the last row
            temp_j = minus_inf;
            for(int l = 1; l < L2; l++)
            {
                if( F[L1-1][l] > temp_j )
                {
                    j_max = l;
                    temp_j = F[L1-1][l];
                }
            }

            if (temp_j > temp_i) // Start with the max of the last row
            {
              i_max = L1-1;
              while(j>j_max)
              {
                  seq_1_al = '-' + seq_1_al;
                  seq_2_al = seq_2[j-1] + seq_2_al;
                  j--;
              }
            }
            else                // Start with the max of the last column
            {
              j_max = L2-1;
              while(i>i_max)
              {
                  seq_1_al = seq_1[i-1] + seq_1_al;
                  seq_2_al = '-' + seq_2_al;
                  i--;
              }
            }
        }else{ // No free trailing
            i_max = L1-1;
            j_max = L2-1;
        }

        score = F[i_max][j_max];


        // TRANSVERSAL TRACEBACK
        while( i > 0 && j > 0 )
        {
            if(i > 0 && j > 0 && F[i][j] == F[i-1][j-1] + s[m[seq_1[i-1]]][m[seq_2[j-1]]])
            {
                seq_1_al = seq_1[i-1] + seq_1_al;
                seq_2_al = seq_2[j-1] + seq_2_al;
                i--;
                j--;
            }
            else if(i > 0 && F[i][j] == F[ i-1 ][ j ] - gap_penalty )
            {
                seq_1_al = seq_1[i-1] + seq_1_al;
                seq_2_al = '-' + seq_2_al;
                i--;
            }
            else
            {
                seq_1_al = '-' + seq_1_al;
                seq_2_al = seq_2[j-1] + seq_2_al;
                j--;
            }
        }

        if(j==0)
        {
            while(i>0)
            {
                seq_1_al = seq_1[i-1] + seq_1_al;
                seq_2_al = '-' + seq_2_al;
                i--;
            }
        }
        else
        {
            while(j>0)
            {
                seq_1_al = '-' + seq_1_al;
                seq_2_al = seq_2[j-1] + seq_2_al;
                j--;
            }
        }

        // formating of the score
        string Score;
        ostringstream convert;
        convert << score;
        Score = convert.str();


        string alignment = seq_1_al + '|' + seq_2_al + '|' + Score ;
        //De-allocate memory if terminate
        if (terminate == 1)
        {
        for (int i = 0; i < L1; ++i)
            delete [] F[i];
        delete [] F;
        }

        return alignment;
}
