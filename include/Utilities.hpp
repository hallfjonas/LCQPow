// TODO: Write header


#ifndef LCQPOASES_UTILITIES_HPP
#define LCQPOASES_UTILITIES_HPP

namespace lcqpOASES {
    
    enum returnValue {
        SUCCESSFUL_RETURN = 0,						/**< Successful return. */
        ILLEGAL_ARGUMENT = 1,                       /**< Illegal arguments. */
        LCQPOBJECT_NOT_SETUP = 2,                   /**< Constructor has not been called. */
        INDEX_OUT_OF_BOUNDS = 3,                    /**< Index out of bounds. */
        SUBPROBLEM_SOLVER_ERROR = 4,
        UNABLE_TO_READ_FILE = 5,
        INVALID_PENALTY_UPDATE_VALUE = 6,
        INVALID_COMPLEMENTARITY_TOLERANCE = 7,
        INVALID_INITIAL_PENALTY_VALUE = 8,
        NOT_YET_IMPLEMENTED = 1000                 /**< Not yet implemented (internal use only). */
    };

    enum printLevel {
        NONE = 0,                                   /**< No Output. */                    
        OUTER_LOOP_ITERATES = 1,                    /**< Print stats for each outer loop iterate. */
        INNER_LOOP_ITERATES = 2,                    /**< Print stats for each inner loop iterate. */
        VERBOSE = 3                                 /**< Print stats for each inner loop (and possibly output of subproblem solver). */
    };

    class Options {

        public:
            /** Default constructor. */
            Options( );

            /** Copy constructor (deep copy). */
            Options(	const Options& rhs			/**< Rhs object. */
                        );

            /** Destructor. */
            ~Options( );

            /** Assignment operator. */
            Options& operator=( const Options& rhs );

            void setToDefault( );                   /**< Sets all options to default values. */

            returnValue ensureConsistency( );       /**< Ensures the consistency of given options. */

            double complementarityTolerance;		/**< Complementarity tolerance. */
            double initialComplementarityPenalty;	/**< Start value for complementarity penalty term. */
            double complementarityPenaltyUpdate;	/**< Factor for updating penaltised complementarity term. */

            bool solveZeroPenaltyFirst;             /**< Flag indicating whether first QP should ignore penalization. */

            printLevel printLvl;                    /**< Print level. */

        protected:
            void copy( const Options& rhs );        /**< Copy each property. */
    };


    class Utilities {
        public:
            // C = A*B.
            static void MatrixMultiplication(const double* const A, const double* const B, double* C, int m, int n, int p);

            // C = A'*B + B'*A
            static void MatrixSymmetrizationProduct(const double* const A, const double* const B, double* C, int m, int n);

            // Read integral data from file
            static returnValue readFromFile(int* data, int n, const char* datafilename);

            // Read float data from file
            static returnValue readFromFile(double* data, int n, const char* datafilename );

            constexpr static const double EPS = 1.11e-16;
    };
}


#endif  // LCQPOASES_UTILITIES_HPP