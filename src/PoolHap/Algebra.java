package PoolHap;

import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;

public class Algebra {
    public static double[] minus(double[] a1, double[] a2) {
        if (a1.length != a2.length) {
            System.out.println("Can't carry out minus: a1.length!=a2.length");
            return null;
        }
        double[] result = new double[a1.length];
        for (int k = 0; k < a1.length; k++) {
            result[k] = a1[k] - a2[k];
        }
        return result;
    }

    public static double[] minus(double[] a1, int[] a2) {
        if (a1.length != a2.length) {
            System.out.println("Can't carry out minus: a1.length!=a2.length");
            return null;
        }
        double[] result = new double[a1.length];
        for (int k = 0; k < a1.length; k++) {
            result[k] = a1[k] - a2[k];
        }
        return result;
    }

    public static double[] minus(int[] a1, double[] a2) {
        if (a1.length != a2.length) {
            System.out.println("Can't carry out minus: a1.length!=a2.length");
            return null;
        }
        double[] result = new double[a1.length];
        for (int k = 0; k < a1.length; k++) {
            result[k] = a1[k] - a2[k];
        }
        return result;
    }

    public static String intarray2string(int[] array) {
        String the_str = "";
        for (int k = 0; k < array.length; k++) {
            the_str = the_str + array[k];
        }
        return the_str;
    }

    public static double sum(double[] a) {
//    	System.out.println("\n___________");
        double result = 0;
        for (int k = 0; k < a.length; k++) {
//        	System.out.println(a[k]+"\t"+result);
            result += a[k];
        }
        return result;
    }

    public static double sum(ArrayList<Double> a) {
        double result = 0;
        for (int k = 0; k < a.size(); k++) {
            result += a.get(k);
        }
        return result;
    }

    public static double[] add(double[] a1, double[] a2) {
        if (a1.length != a2.length) {
            System.out.println("Can't carry out add: a1.length!=a2.length");
            return null;
        }
        double[] result = new double[a1.length];
        for (int k = 0; k < a1.length; k++) {
            result[k] = a1[k] + a2[k];
        }
        return result;
    }

    public static double[] add(double[] a1, double c) {
        double[] result = new double[a1.length];
        for (int k = 0; k < a1.length; k++) {
            result[k] = a1[k] + c;
        }
        return result;
    }

    public static double[] times(double[] a1, double[] a2) {
        if (a1.length != a2.length) {
            System.out.println("Can't carry out inner mutiplication: a1.length!=a2.length");
            return null;
        }
        double[] result = new double[a1.length];
        for (int k = 0; k < a1.length; k++) {
            result[k] = a1[k] * a2[k];
        }
        return result;
    }

    public static double[] times(double[] a1, double c) {
        double[] result = new double[a1.length];
        for (int k = 0; k < a1.length; k++) {
            result[k] = a1[k] * c;
        }
        return result;
    }

    public static double[][] times(double[][] a1, double c) {
        double[][] result = new double[a1.length][];
        for (int k = 0; k < a1.length; k++) {
            result[k] = times(a1[k], c);
        }
        return result;
    }

    public static double[] divide(double[] a1, double[] a2) {
        if (a1.length != a2.length) {
            System.out.println("Can't carry out inner divide: a1.length!=a2.length");
            return null;
        }
        double[] result = new double[a1.length];
        for (int k = 0; k < a1.length; k++) {
            result[k] = a1[k] / a2[k];
        }
        return result;
    }

    public static double inner_product(double[] a1, double[] a2) {
        if (a1.length != a2.length) {
            System.out.println("Can't carry out inner product: a1.length!=a2.length");
            return Double.NaN;
        }
        double result = 0;
        for (int k = 0; k < a1.length; k++) {
            result += a1[k] * a2[k];
        }
        return result;
    }

    public static double mean(double[] p) {
        return sum(p) / p.length;
    }

    public static void rmlow_and_normalize(double[] p, double zero_cutoff) {
        for (int k = 0; k < p.length; k++) {
            if (p[k] < zero_cutoff) {
                p[k] = 0;
            }
        }
        double sum = sum(p);
        for (int k = 0; k < p.length; k++) {
            p[k] = p[k] / sum;
        }
    }

    /*
     *  Given an array, divide all the entries by their sum to form a distribution
     */
    public static void normalize_ditribution(double[] p) {

    	
        double sum = sum(p);
        for (int k = 0; k < p.length; k++) {
            p[k] = p[k] / sum;
        }

    }

    public static void normalize_ditribution(ArrayList<Double> p) {
        double sum = sum(p);
        for (int k = 0; k < p.size(); k++) {
            double x = p.get(k) / sum;
            p.set(k, x);
        }
    }

    /*
     * assigning NaN with an average. the overall sum of all NaN entries will the the "proportion"
     * of others
     */
    public static void normalize_ditribution(ArrayList<Double> p, double proportion) {
        double sum0 = 0;
        int total_NaN = 0;
        for (int k = 0; k < p.size(); k++) {
            if (!Double.isNaN(p.get(k))) {
                sum0 += p.get(k);
            }
            else total_NaN++;
        }
        if (sum0 == 0) { // all are NaN
            for (int k = 0; k < p.size(); k++) {
                p.set(k, 1.0 / p.size());
            }
            return;
        }
        double sum = sum0 / (1 - proportion);
        double nan_average = sum * proportion / total_NaN;
        for (int k = 0; k < p.size(); k++) {
            if (!Double.isNaN(p.get(k))) {
                double x = p.get(k) / sum;
                p.set(k, x);
            } else {
                p.set(k, nan_average);
            }
        }
    }

    public static String intarray2string(int[] a, String sep) {
        String st = "";
        for (int k = 0; k < a.length - 1; k++) {
            st = st + a[k] + sep;
        }
        st = st + a[a.length - 1];
        return st;
    }

    public static double[] abs(double[] p) {
        double[] abs = new double[p.length];
        for (int k = 0; k < p.length; k++) {
            abs[k] = Math.abs(p[k]);
        }
        return abs;
    }

    public static double[] exp(double[] a) {
        double[] exp = new double[a.length];
        for (int k = 0; k < a.length; k++) {
            exp[k] = Math.exp(a[k]);
        }
        return exp;
    }

    public static double quadratic_form(double[] X, double[][] A) {
        RealMatrix AA = MatrixUtils.createRealMatrix(A);
        double[] tXA = AA.preMultiply(X);
        return inner_product(tXA,X);
    }

    public static double[] diag(double[][] A) {
        double[] diag = new double[A.length];
        for (int k = 0; k < A.length; k++) {
            diag[k] = A[k][k];
        }
        return diag;
    }

    /*
     *  density of normal distribution as likelihood
     */
    public static double logL_normal(double[][] sigma, double[] mu, double[] x) {
        int n = mu.length;
        if (x.length != n || sigma.length != n || sigma[0].length != n) {
            System.out.println(x);
        }
        double log_likelihood = Math.log(2 * Math.PI) * (-n / 2.0);
        SingularValueDecomposition svd = new SingularValueDecomposition(
            MatrixUtils.createRealMatrix(sigma));

        double[] sv = svd.getSingularValues();
        double log_p_det = 0;
        for (int k = 0; k < sv.length; k++) {
            if (sv[k] != 0) {
                log_p_det = log_p_det + Math.log(sv[k]);
            }
        }
        log_likelihood = log_likelihood - 0.5 * log_p_det;
        RealMatrix si = svd.getSolver().getInverse();
        double quadratic = quadratic_form(minus(x, mu), si.getData());
        log_likelihood = log_likelihood + (-0.5 * quadratic);
        return log_likelihood;
    }

    public static double distance(double[] x, double[] y) {
        double dis = 0;
        for (int k = 0; k < x.length; k++) {
            dis += Math.abs(x[k] - y[k]);
        }
        return dis;
    }

    public static boolean not_in_span(
        SingularValueDecomposition svd,
        double[][] sigma,
        double[] mu,
        double[][] x) {

        double minimal = 0.001;
        DecompositionSolver solver = svd.getSolver();
        for (int k = 0; k < x.length; k++) {
            double[] B = minus(mu, x[k]);
            double[] X = solver.solve(new ArrayRealVector(B)).toArray();
            double[] Bprime = (new Array2DRowRealMatrix(sigma)).operate(X);
            double dis = distance(B, Bprime);
            if (dis > minimal) {
                return true;
            }
        }
        return false;
    }

    /*
     *  product of multiple density sharing the same normal distribution (\sigma and \mu)
     *  To save time of solving the matrix we can't run the above method k times.
     */
    public static double logL_aems(double[][] sigma, double[] mu, double[][] x) {
        // Formerly, logL_normal.
        // x is public double[][] inpool_site_freqs i.e.: #loci x #pools, so it should be x.length
        int n = mu.length;
        int dim_span = 0;

        // TODO: [LEFTOVER]
        // System.out.println(x.length + "\t" + sigma.length + "\t" + sigma[0].length);

        if (x.length != n || sigma.length != n || sigma[0].length != n) {

            // If the number of variant positions don't match up between all three input parameters,
            // report the error.
            System.out.println("logL_aems: Length of arrays don't match!");
        }
        SingularValueDecomposition svd = new SingularValueDecomposition(
            MatrixUtils.createRealMatrix(sigma));

        // Solve the covariate matrix using singular value decomposition to get its eigenvalues.
        double[] sv = svd.getSingularValues();

        // TODO: [LEFTOVER]
        // if (not_in_span(svd, sigma, mu, x)) {
        //     return Double.NaN; // (0-i) If the the haplotypes are not in the span return an NaN.
        // }

        double log_p_det = 0;
        for (int k = 0; k < sv.length; k++) {
            if (sv[k] != 0) {
                // If the singular value/eigenvalue at variant k is not 0, add log(SV) to log_p_det
                // (starts at 0). Also, increment dim_span.
                log_p_det = log_p_det + Math.log(sv[k]);
                dim_span++;
            }
        }

        // TODO: [LEFTOVER]
        // System.out.println("dim_span: " + dim_span);

        // Instead of using the number of variants (mu.length), Set the base value of the
        // log-likelihood (LL) as...
        double constant = Math.log(2 * Math.PI) * (-dim_span / 2.0);

        // ... log(2 * pi) / (1/2 * size_min_set_span), and subtract half of log_p_det from the
        // log-likelihood.Recall that log_p_det is the sum of the log(covariate eigenvalues).
        constant = constant - 0.5 * log_p_det;

        // Get the inverse matrix of the covariate eigenvalues.
        RealMatrix si = svd.getSolver().getInverse();
        int num_x = x[0].length;
        double logL = 0;
        for (int k = 0; k < num_x; k++) { // For each pool in the double[][] data matrix...
            double[] data_inpool = new double[n];
            for (int v = 0; v < n; v++) {
                data_inpool[v] = x[v][k];
            }

            // (0-ii) Get the inner product of (data - mu) and (inverse matrix of covariate
            // eigenvalues). Add that base value to the running total of the log-likelihood.
            double quadratic = quadratic_form(minus(data_inpool, mu), si.getData());

            // ...subtract half of the quadratic form from the base value of the log-likelihood.
            double log_likelihood = constant + (-0.5 * quadratic);

            // In other words, the base log-likelihood decreases in a different way for each pool,
            // and total log-likelihood is a sum of all of the individual ones.
            logL = logL + log_likelihood;
        }

        // TODO: [LEFTOVER]
        // System.out.println("debug: " + logL);

        return logL;
    }

    public static double logL_rjmcmc(double[][] sigma, double[] mu, double[] x) {
        // NOTE: This might be wayyyyyy too small to use. Will have to observe the values coming out
        // of it.
        int n = mu.length;
        int dim_span = 0;

        // Note that x is the DATA i.e.: the observed allele frequencies in pool p, as set by
        // update_sigma_mu_logL.
        // If the number of variant positions don't match up between all three input parameters,
        // report the error.
        if (x.length != n || sigma.length != n || sigma[0].length != n){
            System.out.println("logL_rjmcmc: Length of arrays don't match!");
        }
        SingularValueDecomposition svd = new SingularValueDecomposition(
            MatrixUtils.createRealMatrix(sigma));

        // Solve the covariate matrix using singular value decomposition to get its eigenvalues.
        double[] sv = svd.getSingularValues();

        // TODO: [LEFTOVER]
        // if (not_in_span(svd, sigma, mu, x)) {
        //    return Double.NaN; // (0-i) If the the haplotypes are not in the span return an NaN.
        // }

        double log_p_det = 0;
        for (int k = 0; k < sv.length; k++) {
            // If the singular value/eigenvalue at variant k is not 0, add log(SV) to log_p_det
            // (starts at 0). Also, increment dim_span.
            if (sv[k] != 0){
                log_p_det = log_p_det + Math.log(sv[k]);
                dim_span++;
            }
        }

        // TODO: [LEFTOVER]
        //System.out.println("dim_span: " + dim_span);

        // Instead of using the number of variants (mu.length), Set the base value of the
        // log-likelihood (LL) as...
        double constant = Math.log(2 * Math.PI) * (-dim_span / 2.0);

        // ... log(2 * pi) / (1/2 * size_min_set_span), and subtract half of log_p_det from the
        // log-likelihood. Recall that log_p_det is the sum of the log(covariate eigenvalues).
        constant = constant - 0.5 * log_p_det;

        // Get the inverse matrix of the covariate eigenvalues.
        RealMatrix si = svd.getSolver().getInverse();
        double logL = 0;

        // (0-ii) Get the inner product of (data - mu) and (inverse matrix of covariate
        // eigenvalues). Add that base value to the running total of the log-likelihood.
        double quadratic = quadratic_form(minus(x, mu), si.getData());

        // ...subtract half of the quadratic form from the base value of the log-likelihood.
        double log_likelihood = constant + (-0.5 * quadratic);

        // In other words, the base log-likelihood decreases in a different way for each pool, and
        // total log-likelihood is a sum of all of the individual ones.
        logL = logL + log_likelihood;

        // TODO: [LEFTOVER]
        //System.out.println("debug: " + logL);

        return logL;
    }
    
    /**
     * Sort the array double[] frequency, and adjust their original indexes accordingly; 
     * and then, select the top ones and return a boolean array labeling the indexes that 
     * will be removed as TRUE.
     *  
     * Before sorting, random permute the arrays so that order will be disrupted: this is to 
     * introduce some randomness so that, if one take the top ones after sorting,  one will 
     * randomly get different ones when there are many entries with the same values.  
     * 
     * This function will be used to filter low-frequnt haplotypes in AEM and regional LASSO
     * 
     * @param frequency
     * @param num_remained_tops
     */
    public static boolean[] permute_sort_and_remove(double[] frequency, int num_remained_tops) {
        int array_len=frequency.length;
        int[] indexes=new int[array_len];
        for(int k=0;k<array_len;k++) {
            indexes[k]=k;
        }
        // permute n*n times first so that the equal entries are disrupted.
        for(int round=0; round<array_len; round++) {
            for(int i=0; i<array_len; i++) {
                int j = (int) (Math.random()*(array_len-1));
                double tmp=frequency[i];
                frequency[i]=frequency[j];
                frequency[j]=tmp;
                int tmp_index=indexes[i];
                indexes[i]=indexes[j];
                indexes[j]=tmp_index;
            }
        }
        // sort the array and its indexes accordingly:
        for(int i=0; i<array_len-1; i++) {
            for(int j=array_len-1; j>i; j--) {
                if(frequency[j-1]<frequency[j]) { // swap them and the indexes
                    double tmp=frequency[j-1];
                    frequency[j-1]=frequency[j];
                    frequency[j]=tmp;
                    int tmp_index=indexes[j-1];
                    indexes[j-1]=indexes[j];
                    indexes[j]=tmp_index;
                }
            }
        }
//        for(int i=0; i<array_len; i++) {
//            System.out.println(frequency[i]+":\t"+indexes[i]);
//        }
        boolean[] removed_indexes = new boolean[array_len];
        for(int remained=0;remained<num_remained_tops;remained++) {
            removed_indexes[indexes[remained]]=false;
        }for(int removed=num_remained_tops; removed<array_len;removed++) {
            removed_indexes[indexes[removed]]=true;
        }
//        for(int i=0; i<array_len; i++) {
//            if(removed_indexes[i])
//                System.out.print(i+" ");
//        }
        return removed_indexes;
    }
    
//    public static void main(String[] args) {
//        double[] array= {0.9,0.1,0.3,0.3,0.1,0.1,0.7,0.8,0.4,0.9,0.6,0.5,0.9,0.1};
//        int remained_tops=9;
//        permute_and_sort(array, remained_tops);
//    }
}
