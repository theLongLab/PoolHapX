package MiscFunctions;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Map;
import java.util.ArrayList;
import java.util.HashMap;

public class PairedReadLinker {
	
//	public static void link_paired_vef(String input_raw_vef, String output_vef) throws IOException {
//		HashMap<String,Integer> dict_read = new HashMap<String,Integer>();
//		ArrayList<String>  readname_arr= new ArrayList<String >();
//		ArrayList<String>  pos_arr= new ArrayList<String >();
//		ArrayList<Integer>  readstart_arr= new ArrayList<Integer >();
//		ArrayList<Integer>  readend_arr= new ArrayList<Integer >();
//		ArrayList<String>  readrange_arr= new ArrayList<String >();
//		BufferedReader bufferedreader= new BufferedReader(new FileReader(input_raw_vef));
//		String line= null;
//		int index=0;
//		
//		while ( (line =bufferedreader.readLine())!=null ){
//			line= line.replace("\r","");
//			String[] line_arr = line.split("\t");
//			String readname = line_arr[0];
//			int readstart = Integer.parseInt(line_arr[3]);
//			
//			
//			if (dict_read.containsKey(readname)) {
//				String[] vars = line_arr[1].split(";|=");
//				String to_add = ""; 
//				for (int v = 0; v < vars.length; v += 2)
//					if (!pos_arr.get(dict_read.get(readname)).contains(vars[v])) {
//						to_add = to_add + vars[v] + "=" + vars[v + 1] + ";"; 
//					}
//				if 	(readstart>=readstart_arr.get(dict_read.get(readname))){
//					readrange_arr.set(dict_read.get(readname), 
//							readrange_arr.get(dict_read.get(readname))+"\t"+line_arr[3]+ "\t"+line_arr[4]); 
//					pos_arr.set(dict_read.get(readname), pos_arr.get(dict_read.get(readname))+to_add);
//				}else {
//					readrange_arr.set(dict_read.get(readname),line_arr[3]+ "\t"+line_arr[4]+"\t" 
//									+readrange_arr.get(dict_read.get(readname))); 
//					pos_arr.set(dict_read.get(readname), to_add+pos_arr.get(dict_read.get(readname)));		
//				}
//				
//			}else {
//				dict_read.put(readname, index);
//				index ++;
//				readname_arr.add(line_arr[0]);
//				pos_arr.add(line_arr[1]);
//				readstart_arr.add(Integer.parseInt(line_arr[3]));
//				readend_arr.add(Integer.parseInt(line_arr[4]));
//				readrange_arr.add(line_arr[3]+"\t"+line_arr[4]);
//			}
//		}
//		bufferedreader.close();
//		
//		PrintWriter pw = new PrintWriter(new FileWriter(output_vef));		
//		for (int i =0;i<readname_arr.size();i++ ) {
//			String tmp_str = readname_arr.get(i)+"\t"+pos_arr.get(i)+"\t"+"//"+"\t"+readrange_arr.get(i)+"\n";
//			pw.write(tmp_str);
//		}
//        pw.flush();
//		pw.close();
//		System.out.println(output_vef+" based on "+input_raw_vef+" has been generated.");
//	}
	
	public static String pos_sort  (String x   )
			throws IOException {
		ArrayList<Integer >  pos_arr = new ArrayList<Integer>();
		ArrayList<String >  allele_arr = new ArrayList<String>();
		String[] x_arr = x.split(";");
		
		for (int v = 0; v < x_arr.length; v ++) {
			String tmp = x_arr[v];
			String[] tmp_arr = tmp.split("=");
			int pos = Integer.parseInt(tmp_arr[0]);
			String allele = tmp_arr[1];
			boolean flag =true;
			for (int i=0; i< pos_arr.size();i++) {
				if (pos_arr.get(i) == pos) {
					flag =false;
					break;
				}
			}
			if (flag) {
				pos_arr.add( pos);
				allele_arr.add(allele); 
			}
		}
		for (int i= 0;i< pos_arr.size();i++) {
			for (int j= i;j< pos_arr.size();j++) {
				if (pos_arr.get(i) > pos_arr.get(j) ) {
					int tmp_int= pos_arr.get(i);
					pos_arr.set(i, pos_arr.get(j));
					pos_arr.set(j,  tmp_int );
					String tmp_str= allele_arr.get(i);
					allele_arr.set(i, allele_arr.get(j));
					allele_arr.set(j,  tmp_str );
				}
			}
		}
		String out="";
		for (int i= 0;i< pos_arr.size();i++) {
			out= out+Integer.toString(pos_arr.get(i))+"="+ allele_arr.get(i)+";";
		}
		return out;
	}
	
	public static void link_paired_vef(String input_raw_vef, String output_vef) throws IOException {
		BufferedReader bufferedreader= new BufferedReader(new FileReader(input_raw_vef));
		String line= null;
		
		HashMap<String,String> read_dict = new HashMap<String,String>();
		ArrayList<String >  read_name = new ArrayList<String>();
		
		
		while ( (line =bufferedreader.readLine())!=null ){
			line= line.replace("\r","");
			String[] line_arr = line.split("\t");
			String readname = line_arr[0];
			String allele= line_arr[1];
			if (read_dict.containsKey(readname)) {
				read_dict.put(readname ,read_dict.get(readname)+  allele ); 
			} else {
				read_dict.put(readname , allele ); 
				read_name.add(readname);
			}
		}
		bufferedreader.close();
		
		PrintWriter pw = new PrintWriter(new FileWriter(output_vef));	
		for (int i=0; i< read_name.size();i++) {
			String tmp_str = read_name.get(i)+"\t"+pos_sort( read_dict.get(read_name.get(i))  ) +"\n";
			pw.write(tmp_str);
		}
		
        pw.flush();
		pw.close();
		System.out.println(output_vef+" based on "+input_raw_vef+" has been generated.");
	}
	
	
//	public static void main(String[] args) {	
//		String prefix = args[0];
//		
//	}

}
