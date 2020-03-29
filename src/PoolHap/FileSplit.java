package PoolHap;
import java.io.*;
import java.util.*;

public class FileSplit {
	public int num_snp_file;

	public FileSplit(String[] vef_arr, String gs_var_pos,  int num_snp_job 
    		) throws IOException {
		this.num_snp_file =  num_snp_job;
		BufferedReader br = new BufferedReader(new FileReader(gs_var_pos));
        String line = "";
        int line_count =-1;
        String header_gs_var_pos = "";
        ArrayList<String> file_arr= new ArrayList<String>();
        ArrayList<String> pos_arr= new ArrayList<String>();
        ArrayList<ArrayList<Integer>> job_index_arr = new ArrayList<ArrayList<Integer>>();
        ArrayList<HashSet<String>> job_snp_dict = new ArrayList<HashSet<String>>();
        while ((line = br.readLine()) != null) {
        	line = line.replace("\r", ""); // remove newline characters
        	if (line_count==-1) {
        		header_gs_var_pos = line;
        	}else {
        		String[] line_arr = line.split("\t");
        		String tmp= line_arr[0];
        		String[] tmp_arr = tmp.split(";");
        		String pos = tmp_arr[1];
        		file_arr.add(line);
        		pos_arr.add(pos);
        		if ((line_count % this.num_snp_file ) ==0) {
        			job_snp_dict.add( new HashSet<String>() );
        			job_index_arr.add( new ArrayList<Integer>());
        		}
        		job_snp_dict.get(job_snp_dict.size()-1).add(pos);
        		job_index_arr.get(job_index_arr.size()-1).add(line_count);
        	}
        	line_count++ ;
        }
        br.close();
        
        for (int i=0; i< job_index_arr.size();i++ ) {
        	FileWriter mydata = new FileWriter(gs_var_pos+ "_"+ Integer.toString(i),false);
            PrintWriter pw = new PrintWriter(mydata);
            pw.write( header_gs_var_pos + "\n"); 
            for (int j=0; j< job_index_arr.get(i).size();j++ ) {
            	pw.write(file_arr.get(job_index_arr.get(i).get(j))+ "\n" );
            }
            pw.flush();
            pw.close();
        }
        
        VefFileSplit(vef_arr, job_snp_dict);
	}
	
	
	public void VefFileSplit (String[] vef_arr , ArrayList<HashSet<String>> job_snp_dict)
			throws  IOException {
		for (int i=0; i<job_snp_dict.size();i++ ) {
			int lastIndex = vef_arr[0].lastIndexOf("/");
			new File(vef_arr[0].substring(0, lastIndex) + "/" +"../split_vef/"+Integer.toString(i) ).mkdir();
		}
		
		for (int i=0; i<job_snp_dict.size();i++ ) {
			for (int j=0;j< vef_arr.length; j++) {
				int lastIndex = vef_arr[j].lastIndexOf("/");
				String path= vef_arr[j].substring(0, lastIndex) + "/" +"../split_vef"+ "/" + 
						Integer.toString(i) + "/" +vef_arr[j].substring(lastIndex);
				FileWriter mydata = new FileWriter(path ,false);
	            PrintWriter pw = new PrintWriter(mydata);
	            BufferedReader br = new BufferedReader(new FileReader(vef_arr [j]));
	            String line = "";
	            while ((line = br.readLine()) != null) {
	            	line = line.replace("\r", ""); // remove newline characters
	                String[] line_arr = line.split("\t"); // Read_Name\tPos=Allele;...\t//\tStart\tEnd
	                if (line_arr[1].contains("=")) {
	                    String geno = line_arr[1]; // the segregating site information column
	                    String [] loci_arr= geno.split(";");
	                    String write_geno= "";
	                    for (int k =0; k< loci_arr.length; k++) {
	                    	String loci= loci_arr[k].split("=")[0];
	                    	if  (job_snp_dict.get(i).contains(loci)) {
	                    		write_geno = write_geno+ loci_arr[k]+";";
	                    	}
	                    }
	                    if (!write_geno.equals("")){
	                    	String write_line = line_arr[0]+"\t" + write_geno;
	                    	for (int k =2; k< line_arr.length; k++) {
	                    		write_line=write_line+ "\t"+ line_arr[k];
	                    	}
	                    	pw.write(write_line+ "\n" );
	                    }
	                }else {
	                	
	                }
	            }
	            br.close();
	            pw.flush();
	            pw.close();
			}
		}
	}
	
}
