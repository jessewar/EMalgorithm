import java.util.ArrayList;
import java.util.Collections;
import java.util.Scanner;

/**
 * Jesse Warren
 * Section AC
 * jessewar
 * HW9
 * 
 */

// Implementation of the EM algorithm where we are given that tau == 1/3 and var == 1
public class EMalgorithm {

	// executes our EM function
	public static void main(String[] args) {
		expectationMaximization();
	}
	
	public static final double epsilon = 0.001;  // constant for determining when to stop iterating
	
	// Reads in a sequence of numbers from standard in and prints out the value of mu1, mu2, mu3, and the loglihood after each iteration
	// Prints the values of z_ij after the final iteration
	public static void expectationMaximization() {
	       ArrayList<Integer> data = new ArrayList<Integer>();
	       Scanner console = new Scanner(System.in);
	       
	       // add all x_i's to data
	       while (console.hasNext()) {
	    	   	 String next = console.next();
	    	   	 if (!next.equals("q")) { data.add(Integer.parseInt(next)); }
	    	   	 else { break; }
	       }
	       Collections.sort(data);
	       System.out.println("Mu_1" + "\t" + "Mu_2" + "\t" + "Mu_3" + "\t" + "Loglihood");
	       // INITIALIZATION CRITERION: mu1 = smallest datum, mu2 = median datum, mu3 = largest datum
	       expectationMaximizationHelper(data, data.get(0), data.get(data.size() / 2), data.get(data.size() - 1));
	}

	// helper function for expectationMaximization()
	private static void expectationMaximizationHelper(ArrayList<Integer> data, double mu1, double mu2, double mu3) {
	       // create array to store which distribution each z_ij value, for j = 1, 2, or 3, has highest expected value
	       int[] z_ijList = new int[data.size()];

	       // E-STEP

	       double logLihood = 0;
	       // calculate, for each x_i, the corresponding z_ij values and store the distribution for which z_ij is highest
	       for (int i = 0; i < data.size(); i++) {
	    	   int x_i = data.get(i);
	           double z_i1 = bayes(mu1, mu2, mu3, x_i);
	           double z_i2 = bayes(mu2, mu1, mu3, x_i);
	           double z_i3 = bayes(mu3, mu1, mu2, x_i);

	           logLihood += Math.log(likelihood(mu1, mu2, mu3, x_i));
	           double max = Math.max(Math.max(z_i1, z_i2), z_i3);
	           z_ijList[i] = max == z_i1 ? 1 : (max == z_i2 ? 2 : 3); // 1 if z_i1 is max, 2 if z_i2 is max, 3 otherwise
	       }
	       
	       // M-STEP

	       double newMu1 = mStep(data, z_ijList, 1);
	       double newMu2 = mStep(data, z_ijList, 2);
	       double newMu3 = mStep(data, z_ijList, 3);

	       // TERMINATION CRITERION: all means differ from last means by less than epsilon == 0.001
	       // if values converge, just print, otherwise print and continue recursive calls with new mu values
	       if (Math.abs(mu1 - newMu1) < epsilon && Math.abs(mu2 - newMu2) < epsilon && Math.abs(mu3 - newMu3) < epsilon) {
	           System.out.println(newMu1 + "\t" + newMu2 + "\t" + newMu3 + "\t" + logLihood);
	           
	           System.out.println();
	           System.out.println("i" + "\t" + "x_i" + "\t" + "z_i1" + "\t" + "z_i2" + "\t" + "z_i3");
		       for (int i = 0; i < data.size(); i++) {
		    	   int x_i = data.get(i);
		           double z_i1 = bayes(mu1, mu2, mu3, x_i);
		           double z_i2 = bayes(mu2, mu1, mu3, x_i);
		           double z_i3 = bayes(mu3, mu1, mu2, x_i);
		           
		           System.out.println(i + "\t" + x_i + "\t" + z_i1 + "\t" + z_i2 + "\t" + z_i3);
		       }
	       } else {
	    	   System.out.println(newMu1 + "\t" + newMu2 + "\t" + newMu3 + "\t" + logLihood);
	           expectationMaximizationHelper(data, newMu1, newMu2, newMu3);
	       }
	}

	// Uses the equation derived in part 2 of my derivations sheet to estimate mu by finding which mu maximizes the likelihood our x_i's and z_ij's
	private static double mStep(ArrayList<Integer> data, int[] z_ijList, int target) {
	        double weightedSum = 0;
	        int count = 0;
	        for (int i = 0; i < data.size(); i++) {
	            if (z_ijList[i] == target) {	
	               weightedSum += data.get(i);
	               count++;
	            }
	        }

	        return weightedSum / count;	// from my derivation: average, weighted by subpop
	}
	
	// computes the individual likelihood of a single x_i datum given parameters mu1, mu2, and mu3
	private static double likelihood(double mu1, double mu2, double mu3, int x_i) {
	       double p1 = getPDFValue(x_i, mu1, 1);
	       double p2 = getPDFValue(x_i, mu2, 1);
	       double p3 = getPDFValue(x_i, mu3, 1);

	       return (p1 + p2 + p3) / 3;  // divide by 3 b/c tau == 1/3 for each distribution
	}
	
	// computes the value of z_ij using the formula I derived in part 1 of my derivation sheet
	private static double bayes(double mu1_target, double mu2, double mu3, int x_i) {
	       double p1 = getPDFValue(x_i, mu1_target, 1);
	       double p2 = getPDFValue(x_i, mu2, 1);
	       double p3 = getPDFValue(x_i, mu3, 1);

	       return p1 / (p1 + p2 + p3);	// from my derivation
	}

	// computes the f(x_i) for the normal distribution with mean = mu and standard deviation = sd
	private static double getPDFValue(int x_i, double mu, double sd) {
	       return (1 / (sd * Math.sqrt(2 * Math.PI))) * (Math.pow(Math.E, -1 * (Math.pow(x_i - mu, 2) / (2 * Math.pow(sd, 2)))));
	}
}
