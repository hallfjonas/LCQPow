// TODO: Include header

#include "Utilities.hpp"
#include <iostream>

namespace lcqpOASES {
    /*
    *	O p t i o n s
    */
    Options::Options( )
    {
        setToDefault( );
    }


    /*
    *	O p t i o n s
    */
    Options::Options( const Options& rhs)
    {
        copy( rhs );
    }


    /*
    *	~ O p t i o n s
    */
    Options::~Options( )
    { }


    /*
    *	o p e r a t o r =
    */
    Options& Options::operator=( const Options& rhs )
    {
        if ( this != &rhs )
        {
            copy( rhs );
        }

        return *this;
    }


    /*
     *   c o p y
     */
    void Options::copy( const Options& rhs ) {

        complementarityTolerance = rhs.complementarityTolerance;
        initialComplementarityPenalty = rhs.initialComplementarityPenalty;
        complementarityPenaltyUpdate = rhs.complementarityPenaltyUpdate;
        solveZeroPenaltyFirst = rhs.solveZeroPenaltyFirst;
        printLvl = rhs.printLvl;
    }


    /*
     *  e n s u r e C o n s i s t e n c y
     */
    returnValue Options::ensureConsistency( ) {
        
        if (complementarityPenaltyUpdate <= 1)
            throw INVALID_PENALTY_UPDATE_VALUE;

        if (complementarityTolerance < Utilities::EPS)
            throw INVALID_COMPLEMENTARITY_TOLERANCE;

        if (initialComplementarityPenalty <= 0)
            throw INVALID_INITIAL_PENALTY_VALUE;

        return SUCCESSFUL_RETURN;
    }


    /*
     *   s e t T o D e f a u l t
     */
    void Options::setToDefault( ) {
        complementarityTolerance      = 1.0e3 * Utilities::EPS;
        initialComplementarityPenalty = 0.01;
    	complementarityPenaltyUpdate  = 2.0;

        solveZeroPenaltyFirst = true;

        printLvl = printLevel::INNER_LOOP_ITERATES;

    }

    /*
     *   M a t r i x M u l t i p l i c a t i o n
     */
    void Utilities::MatrixMultiplication(const double* const A, const double* const B, double* C, int m, int n, int p) {
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < p; j++) {
                C[i*m + j] = 0;
                for (int k = 0; k < n; k++) {
                    C[i*m + j] += A[i*n + k]*B[k*p + j]; 
                }                
            }
        }
    }


    /*
     *  M a t r i x S y m m e t r i z a t i o n P r o d u c t
     */
    void Utilities::MatrixSymmetrizationProduct(const double* const A, const double* const B, double* C, int m, int n) {
        for (int i = 0; i < m; i++) {
            for (int j = 0; j <= i; j++) {
                C[i*m + j] = 0;
                for (int k = 0; k < n; k++) {
                    C[i*m + j] += A[k*n + i]*B[k*n + j] + B[k*n + i]*A[k*n + j]; 
                }                

                // Make symmetric
                C[j*m + i] += C[i*m + j];
            }
        }
    }

    /*
     *	r e a d F r o m F i l e
     */
    returnValue Utilities::readFromFile( int* data, int n, const char* datafilename )
    {
        int i;
        FILE* datafile;

        /* 1) Open file. */
        if ( ( datafile = fopen( datafilename, "r" ) ) == 0 )
        {
            fclose( datafile );
            return returnValue::UNABLE_TO_READ_FILE;
        }

        /* 2) Read data from file. */
        for( i=0; i<n; ++i )
        {
            if ( fscanf( datafile, "%d\n", &(data[i]) ) == 0 )
            {
                fclose( datafile );
                return returnValue::UNABLE_TO_READ_FILE;
            }
        }

        /* 3) Close file. */
        fclose( datafile );

        return SUCCESSFUL_RETURN;
    }


    /*
     *	r e a d F r o m F i l e
     */
    returnValue Utilities::readFromFile( double* data, int n, const char* datafilename )
    {
        int i;
        FILE* datafile;

        /* 1) Open file. */
        if ( ( datafile = fopen( datafilename, "r" ) ) == 0 )
        {
            fclose( datafile );
            return returnValue::UNABLE_TO_READ_FILE;
        }

        /* 2) Read data from file. */
        for( i=0; i<n; ++i )
        {
            if ( fscanf( datafile, "%lf\n", &(data[i]) ) == 0 )
            {
                fclose( datafile );
                return returnValue::UNABLE_TO_READ_FILE;
            }
        }

        /* 3) Close file. */
        fclose( datafile );

        return SUCCESSFUL_RETURN;
    }
}

