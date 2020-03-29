package Viral_Reconstructions_Tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;

public class Hippo_transfer_output {

	public static void main(String[] args) {
		String results_out="results.out";
		String project_name=args[0];
		String vars_file = project_name + "_vars.intra_freq.txt";
		String final_results_out=project_name + ".inter_freq_vars.txt";
		transfer_file(results_out,vars_file, final_results_out);
	}
	
	/*
	 * Read the results.out file, each line contains the haplotype and the frequence of each haplotye
	 * each_line_arry[0] is the haplotype
	 * each_line_arry[1] is the frequency
	 * When each_line_arry[0] is the haplotype, save those loci that is "1" into loci_reads_listlist,
	 * and save all the loci position that could be "1" to haplo_diff_loci.
	 * When each_line_arry[1], each haplotype frequence for each haplotypes into haplo_frequency_list
	 * Write the final_output_file using the sorted haplo_frequency_list, haplo_diff_loci and loci_reads_listlist
	 */
	public static void transfer_file(String results_out, String vars_file, String final_results_out) {
		try {
			BufferedReader br_file=new BufferedReader(new FileReader(results_out));
			BufferedReader br_vars_file=new BufferedReader(new FileReader(vars_file));
			BufferedWriter bw_file=new BufferedWriter(new FileWriter(final_results_out));
//			ArrayList<Integer> haplo_diff_loci=new ArrayList<Integer>();
			ArrayList<String> haplo_frequency_list= new ArrayList<String>();
			ArrayList<String> loci_list= new ArrayList<String>();
			ArrayList<ArrayList<Integer>> loci_reads_listlist=new ArrayList<ArrayList<Integer>>();
			String vars_line = br_vars_file.readLine();
			vars_line = br_vars_file.readLine(); //read the second line
			while(vars_line!=null) {
				String[] each_line_arry=vars_line.split("\t");
				loci_list.add(each_line_arry[0]);
				vars_line = br_vars_file.readLine();
			}
			br_vars_file.close();
			String each_line=br_file.readLine(); // read the first line
			while(each_line!=null) {
				String[] each_line_arry=each_line.split(" ");
				String[] haplo_array=each_line_arry[0].split("");
				ArrayList<Integer> loci_reads_row=new ArrayList<Integer>();
				for(int i=0;i<haplo_array.length;i++) {
					if(haplo_array[i].equals("1")) {
						loci_reads_row.add(i);
					}
				}
				loci_reads_listlist.add(loci_reads_row);
				haplo_frequency_list.add(each_line_arry[1]);
				each_line=br_file.readLine();
			}
			br_file.close();
			// write the haplotype ID
			bw_file.write("Hap_ID"+"\t");
			for(int i=0;i<loci_reads_listlist.size()-1;i++) {
				bw_file.write("h" + i +"\t");
			}
			bw_file.write("h" + (loci_reads_listlist.size()-1) +"\n");
			// write the frequency of each haplotype
			bw_file.write("freq"+"\t");
			for(int i=0;i<haplo_frequency_list.size()-1;i++) {
				bw_file.write(haplo_frequency_list.get(i)+"\t");
			}
			bw_file.write(haplo_frequency_list.get(haplo_frequency_list.size()-1)+"\n");
			// Write the haplotype position
			// For each haplotype ID, if at this position equals to 1, write "1", else write "0"
			for(int j=0;j<loci_list.size();j++) {
				bw_file.write(loci_list.get(j)+"\t");
				for(int num_hap=0;num_hap<loci_reads_listlist.size()-1;num_hap++) {
					if(loci_reads_listlist.get(num_hap).contains(j)) {
						bw_file.write("1"+"\t");
					}else {
						bw_file.write("0"+"\t");
					}
				}
				if(loci_reads_listlist.get(loci_reads_listlist.size()-1).contains(j)) {
					bw_file.write("1"+"\n");
				}else {
					bw_file.write("0"+"\n");
				}
			}
			bw_file.close();
		}catch(Exception e) {e.printStackTrace();}
		
	}

}
