package Viral_Reconstructions_Tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;

public class Hippo_datafile_transfer {

	public static void main(String[] args) {
		String project_name=args[0];
		String input_datafile = project_name + "_vars.intra_freq.txt";
		String output_datafile= project_name + "_vars.intra_freq_datafile";
		transfer_datafile(input_datafile,output_datafile);
	}
	
	public static void transfer_datafile(String input_datafile, String output_datafile) {
		try {
			BufferedReader br_file=new BufferedReader(new FileReader(input_datafile));
			BufferedWriter bw_file=new BufferedWriter(new FileWriter(output_datafile));
			ArrayList<ArrayList<String>> loci_reads_listlist=new ArrayList<ArrayList<String>>();
			int poolsize=100;
			String line=br_file.readLine();// read the first line, the first line is the pool_ID
			line=br_file.readLine(); //read from the second line 
			while(line!=null) {
				String[] each_position=line.split("\t");
				ArrayList<String> loci_reads_row=new ArrayList<String>();
				for(int i=1;i<each_position.length;i++) {
					loci_reads_row.add(each_position[i]);
				}
				loci_reads_listlist.add(loci_reads_row);
				line=br_file.readLine();
			}
			br_file.close();
			//start to write file
			for(int i=0;i<loci_reads_listlist.get(0).size();i++) { // for each row that we read from the loci_reads_listlist
				double each_loci=0;
				int each_loci_int=0;
				bw_file.write(poolsize+" "); 
				for(int j=0;j<loci_reads_listlist.size()-1;j++) {
					each_loci= Double.parseDouble(loci_reads_listlist.get(j).get(i))*100;
					each_loci_int =(int)each_loci;	
					bw_file.write(each_loci_int+" ");
				}
				each_loci= Double.parseDouble(loci_reads_listlist.get(loci_reads_listlist.size()-1).get(i))*100;
				each_loci_int=(int)each_loci;
				bw_file.write(each_loci_int+"\n");
			}
			bw_file.close();
		}catch(Exception e) {e.printStackTrace();}
		
	}

}
