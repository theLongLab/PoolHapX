package MiscFunctions;


public class SimpleMath {
    public static int sum(int[] a){
        int result=0;
        for(int k=0;k<a.length;k++)result+=a[k];
        return result;
    }
    
    public static double mean(int[] a) {
        return ((double) sum(a)) / a.length;
    }

    public static double stdev(int[] a) {
        double stdev = 0.0; 
        double mean = mean(a);
        for (int i : a) {
            stdev += Math.pow((double) i - mean, 2);
        }
        return Math.sqrt(stdev/a.length);       
    }
    
    public static double nCr(int n, int r){
        int rfact = 1, nfact = 1, nrfact = 1, temp1 = n - r, temp2 = r;
        if (r > n - r) {
            temp1 = r;
            temp2 = n - r;
        }
        for (int i = 1; i <= n; i++) {
            if (i <= temp2) {
                rfact *= i;
                nrfact *= i;
            } else if(i <= temp1) {
                nrfact *= i;
            }
            nfact *= i;
        }
        return nfact / (double) (rfact * nrfact);
    }
}
