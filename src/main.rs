use std::env;
use std::fs;
use std::collections::HashMap;
use serde::{Serialize, Deserialize};
use std::fs::File;
use std::io::{Write, Read, BufReader, BufRead};
use std::cmp::min;
use std::str::Chars;
use std::iter::Rev;
use std::fmt::Debug;
use rust_htslib::bam::{Writer, Format, Header};

//const ALPHABET : [char; 2] = ['a', 'b']; 
const ALPHABET : [char; 4] = ['A', 'C', 'G', 'T'];

fn main() {
	let args : Vec<String> = env::args().collect();
    
    match args[1].as_str(){
    	"index" => {
    		let input = &args[2];
    		let output= &args[3];

    		index(input, output, None);
    	}
    	,
    	"align" => {
    		let input = &args[2];
    		let reads = &args[3];
    		let output = &args[3];


    		align(input, reads, output);
    	},
    	"substring" => {
    		let input = &args[2];
    		let start = &args[3];
    		let end = &args[4];

    		substring(input, start.parse::<usize>().unwrap(), end.parse::<usize>().unwrap());
    	},
    	_ => {println!("Unrecognized command");},
    }
}

//Struct for serialization
#[derive(Serialize, Deserialize, Debug)]
struct FmIndex {
	name : String,
	len : usize,
	first_column : HashMap::<char, usize>,
	last_column : Vec::<char>,
	occ : HashMap<char, Vec<usize>>,
	suffix_arr : Vec<usize>,
	seq : String
}

struct GenRead{
	name: String,
	id : u32,
	seq: String,
}

struct Alignment {
	score : u32,
	CIGAR : String,
	pos: u32
}

//The name and length of the input sequence (which is referred to below as S) : Done
//An efficient representation of the first column of BWM(S) : Done (not efficient tho)
//The last column of BWM(S) : Done
//The occ table (can be sampled if you want, but it need not be) for the characters ‘A’, ‘C’, ‘G’, and ‘T’ Done
//The suffix array (can be sampled if you want, but it need not be) of S Done

fn parse_reads(lines : Vec<String>) -> Vec<GenRead>{
	let mut reads = Vec::<GenRead>::new();

	for i in (0..lines.len()).step_by(2) {
		let line = lines.get(i).unwrap();

		let mut s = line.split('.');
		let name = s.next().unwrap();
		let id = s.next().unwrap().parse::<u32>().unwrap();

		let line = lines.get(i + 1).unwrap();

		reads.push(GenRead {
			name: name.to_string(),
			id: id,
			seq: line.to_string()
		});
	}

	return reads;
}

impl FmIndex {
	fn get_interval(&self, query : &str) -> ((usize, usize), usize){
		//println!("Query: {}", query);

		//Reverse the string so we can iterate it backwards
		let mut remaining_query = query.chars().rev();
		let mut c1 = remaining_query.next().unwrap();

		//Match length will always be at least 1 because you'll always find the first character
		let mut match_len = 1;

		let preceding_rows = {
			//Start sum @ 1 for the temrinating character
			let mut sum = 1;
			for c in ALPHABET.iter() {
				if c != &c1 {
					sum += self.first_column[c];
				}
				else{
					break;
				}
			}

			sum
		};

		let mut bwt_interval : (usize, usize) = (preceding_rows, preceding_rows + self.first_column[&c1]);
		for c2 in remaining_query {
			//println!("Seraching for {} in range : {:?} ", c2, bwt_interval);
			bwt_interval 
				= match self.interval_step(c2, bwt_interval) {
					None => {println!("mismatch"); return (bwt_interval, match_len)},
					Some(interval) => { match_len+= 1; interval }
				};
		}

		return (bwt_interval, match_len);
		//return self.interval_step(c ,(0, self.first_column[&c] - 1), remaining_query, match_len);
	}

	//Returns None on mismatch
	//Returns interval of [a,b) otherwise
	fn interval_step(&self, search_c : char, (bwt_start, bwt_end) : (usize, usize)) -> Option<(usize, usize)> {
		//Retrieve the startiing rank from previous row (won't be out of bounds due to terminating character)
		let start_rank : usize = self.occ[&search_c][bwt_start - 1];
		let end_occ : usize = self.occ[&search_c][bwt_end - 1];

		//If the entries in the occ_table are the same then there is no change within the range and therefore a mismatch
		if start_rank == end_occ {
			return None;
		}

		//Get the rows precceding search_c
		let preceding_rows = {
			//Start sum @ 1 for the temrinating character
			let mut sum = 1;
			for c in ALPHABET.iter() {
				if c != &search_c {
					sum += self.first_column[c];
				}
				else{
					break;
				}
			}

			sum
		};

		return Some((preceding_rows + start_rank, preceding_rows + end_occ));
	}

	fn ref_positions(&self, (start_interval, end_interval) : (usize, usize), seed_end : usize, match_len : usize) -> Vec<usize>{
		let mut alignments : Vec<usize> = Vec::new();

		//Inclusive range
		for i in start_interval..end_interval {
			//The distance from the start of the match to the start of the read
			let dist_to_read_start = seed_end - match_len;

			//Where the read should start
			let pos = self.suffix_arr[i] - dist_to_read_start;
			alignments.push(pos);
		}

		return alignments;
	}

}

fn print_2d_vec<T: Debug>(matrix : &Vec<Vec<T>>){
	for y in matrix {
		for x in y{
			print!("{:?} ", x);
		}
		print!("\n");
	}
}

//semi global alignment trying to fit 'y' into 'x'
fn fitting_alignment(x : &str, y : &str) {
	//Add one for the gap in the start
	let x_len = x.len() + 1;
	let y_len = y.len() + 1;

	const MATCH_SCORE : i32 = 10;
	const GAP_PENALTY : i32 = -7;
	const MISMATCH_PENALTY : i32 = -5;
	
	let mut matrix : Vec<Vec<i32>> = (0..y_len).map(|_| (0..x_len).map(|_| 0).collect()).collect();
	let mut path_next : Vec<Vec<Vec<(usize, usize)>>> = (0..y_len).map(|_| (0..x_len).map(|_| Vec::<(usize, usize)>::new()).collect()).collect();

	let cost_match = |x_pos, y_pos| {
		//Subtract one to account for the gap
		if x.as_bytes()[x_pos - 1] == y.as_bytes()[y_pos - 1] {
			return MATCH_SCORE;
		}
		else {
			return MISMATCH_PENALTY;
		}
	};

	//Populate the first 'x' row to 0
	for i in 0..x_len {
		matrix[0][i] = 0;
	}

	//Populate the first column
	for i in 0..y_len {
		matrix[i][0] = (i as i32) * GAP_PENALTY;
	}

	//Populate the rest of the matrix while keeping track of where the new value originates from
	for y in 1..y_len {
		for x in 1..x_len {
			//All the possbiel values for the next value from the surrounding 3 cells
			//(value, (y, x))
			let options = vec![
				(matrix[y - 1][x - 1] + cost_match(x, y), (y - 1, x - 1)),
				//Add gap penalty to the vertical and horizontal traversal
				(matrix[y - 1][x] + GAP_PENALTY, (y - 1, x)), 
				(matrix[y][x - 1] + GAP_PENALTY, (y, x - 1))
			];
			
			let mut possible_paths = Vec::<(usize, usize)>::new();
			let mut max_score : i32 = i32::MIN;
			//Find the min of the three possible options
			for i in 0..3 {
				if options[i].0 > max_score {
					max_score = options[i].0;
					possible_paths = vec![options[i].1];
				}
				else if options[i].0 == max_score {
					possible_paths.push(options[i].1);
				}
			}

			//Set the updated value and next path step
			matrix[y][x] = max_score;
			path_next[y][x] = possible_paths;
		}

	}
	print_2d_vec(&matrix);

	println!("");
	print_2d_vec(&path_next);

	//Find the largest value in the last row 
	let mut max_val = i32::MIN;
	//(y,x)
	let mut pos = (0,0);

	for i in 0..x_len {
		//TODO: account for multiple max_vals
		if matrix[y_len - 1][i] > max_val {
			max_val = matrix[y_len - 1][i];
			pos = (y_len - 1, i);
		}
		else if matrix[y_len - 1][i] == max_val{

		}
	}

	println!("Max val {:?} @ {:?}", max_val, pos);

	//Begin the backtracing process
	let paths = &path_next[pos.0][pos.1];

	//pos : (y,x)
	fn get_cigars (paths : &Vec<(usize, usize)>, (pos_y, pos_x) : (usize, usize), cigar_acc : String, x : &str, y : &str, matrix : &Vec<Vec<i32>>, path_next : &Vec<Vec<Vec<(usize, usize)>>>) -> Vec<String>{
		//Base case for when we reach the boundries of the matrix
		if paths.len() == 0 {
			return vec![cigar_acc];
		}

		let updated_cigar = |(curr_y,curr_x) : (usize, usize), (next_y, next_x) : (usize, usize), curr_cigar : &str| -> String {
			//pos_ - paths[0]. will return the direction of the path
			match (curr_y - next_y, curr_x - next_x){
				(1, 1) => {
					if x.as_bytes()[curr_x - 1] == y.as_bytes()[curr_y - 1] {
						curr_cigar.to_string() + "M"
					}
					else {
						curr_cigar.to_string() + "X"
					}
				},
				(1, 0) => {
					curr_cigar.to_string() + "I"
				},
				(0, 1) => {
					curr_cigar.to_string() + "D"
				},
				(_, _) => {
					panic!("Bruh wtf.");
				}
			}
		};

		let mut next = (paths[0].0, paths[0].1);
		//Parse the first possible path
		let mut ret_vec = get_cigars(&path_next[next.0][next.1], next, updated_cigar((pos_y, pos_x), next, &cigar_acc), x, y, &matrix, &path_next);

		for i in 1..paths.len() {
			next = (paths[i].0, paths[i].1);
			ret_vec.append(&mut get_cigars(&path_next[next.0][next.1], next, updated_cigar((pos_y, pos_x), next, &cigar_acc), x, y, &matrix, &path_next));
		}

		return ret_vec;
	};

	println!("{:?}", get_cigars(paths, pos, "".to_string(), x, y, &matrix, &path_next));
}

fn align(bwt_file : &str, reads : &str, output : &str){
	//Read the result of the index command in
	let mut file = File::open(bwt_file).unwrap();
	let mut data = Vec::<u8>::new();
	file.read_to_end(&mut data);
	let fm_index : FmIndex = bincode::deserialize(&data[..]).unwrap();

    let file = File::open(reads).expect("no such file");
    let buf = BufReader::new(file);

    let reads = buf.lines()
        .map(|l| l.expect("Could not parse line"))
        .collect();

	let reads = parse_reads(reads);
	const gap : usize = 5;

	for read in reads.iter() {
		let alignments : Vec<Alignment> = Vec::new();
		let read_len = read.seq.len();
		let mut best_score : u32 = u32::MAX;
		let seed_pos = 0;

		let seed_skip = |len : usize| -> usize {(len as f64 / 5.0).floor() as usize};
		let skip = seed_skip(read_len);
		println!("Read: {:?}", &read.seq);

		for seed_start in (0..read_len).step_by(skip){
			let seed_end = min(read_len, seed_start + skip);
			let seed = &read.seq[seed_start..seed_end];

			//TODO: this is very slow. maybe replace
			if seed.to_string().contains("N") {
				continue;
			}

			println!("Seed: {:?}", seed);

			let (interval, match_len) = fm_index.get_interval(seed);

			for pos in fm_index.ref_positions(interval, seed_end, match_len){
				let pos = pos as usize;
				let read_len = read_len as usize;
				let start = (gap + pos);
				let end = (pos + read_len + gap);

				println!("Aligning {} to {}", read.seq, &fm_index.seq[start..end]);
				//let alignment = fm_index.fitting_alignment(read.seq, pos, gap);
				//Looking for the minimum here, not the max
/*				if alignment.score < best_score {
					best_score = alignment.score;
					alignments = vec![alignment; 1];
				}
				else if alignment.score == best_score {
					alignments.push(alignment);
				}*/
			}




		}
	}

}

#[test]
fn test_fitting_alignment(){	
	fitting_alignment("PANCAKE", "ANAK");
}

#[test]
fn test_write_sam(){	
	Writer::from_path("test.sam", &Header::new(), Format::SAM)l
}


/*#[test]
fn test_get_interval(){
	println!("this is a test");
	//let index = index("", "", Some("ATCGATCG$"));
/*	println!("{:?}", index);
*/	//println!("{:?}", index.get_interval("ATC"));
	//println!("{:?}", index.get_interval("CTA"));

	println!("this is a test");
	let index = index("", "", Some("ATCGATCGGCTATAGA$"));
	println!("{:?}", index);
	println!("{:?}", index.get_interval("AT"));

/*


	let index = index("", "", Some("ACAACT$"));
	println!("{:?}", index);
	//Read: ACAT 
	let (interval, match_len) = index.get_interval("CAT");
	//2 because its the value for seed_end is exclusive
	println!("{:?}", interval);
	println!("{:?}", index.ref_positions(interval, 3, match_len));*/

}*/

/*#[test]
fn test_from_file(){
	println!("Testing from file");
	//let index = index("", "", Some("ATCGATCG$"));
/*	println!("{:?}", index);
*/	//println!("{:?}", index.get_interval("ATC"));
	//println!("{:?}", index.get_interval("CTA"));

	println!("this is a test");
	println!("ATCGATCG$");
	let index = index("", "", Some("ATCGATCGGCTATAGA$"));
	println!("{:?}", index);
	println!("{:?}", index.get_interval("AT"));

	let index = index("", "", Some("ACAACT$"));
	println!("{:?}", index);
	//Read: CAC
	let (interval, match_len) = index.get_interval("CA");
	//2 because its the value for seed_end is exclusive
	println!("{:?}", index.ref_positions(interval, 2, match_len));

}*/

fn index(input_file : &str, output_file : &str, test_seq : Option<&str>) -> FmIndex{
	let mut header : String = String::new();
	let mut sequence = 
			{
				match test_seq {
					None => {
						println!("Reading from \'{}\'. Outputing to \'{}\'", input_file, output_file);

						let mut sequence = fs::read_to_string(input_file)
							.expect("Something went wrong reading the file");

						//Read in the sequence and header data from file
						let mut split = sequence.split("\n");

						let mut sequence : String = String::new(); 

						for (i, item) in split.enumerate() {
							if i == 0 { header.push_str(item) }
							else{
								sequence.push_str(item);
							}
						}
						sequence = sequence.trim().to_string();
						//Add the terminating character after processing
						sequence.push_str("$");

						println!("{}", header);
						sequence
					},
					Some(seq) => seq.to_string()
				}
			};

	//TODO: test code
	//let mut sequence : String = String::from("abaaba");

	let len : usize = sequence.chars().count();

	//Constructing the suffix array
	let mut suffix_arr : Vec<(String, usize)> = vec![("".to_string(), 0); len];


	//Populate the suffix array
	for i in 0..len{
		let mut transform = sequence[i..len].to_string();
		transform.push_str(&sequence[0..i]);

		suffix_arr[i] = (transform, i);
	}

	//Sort the suffix array
	suffix_arr.sort_by(|a, b| a.0.cmp(&b.0) );

	//Store the L and F columns
	let mut L_col = Vec::<char>::new();
	let mut F_col = Vec::<char>::new();
	let mut suffix_arr_off = Vec::<usize>::new();

	for (suffix, index) in suffix_arr {
		//TODO: idk for some reason .char_at() is unsafe so i guess i'll just use vectors from now on
		F_col.push(suffix.as_bytes()[0] as char);
		L_col.push(suffix.as_bytes()[len - 1] as char);
		suffix_arr_off.push(index);
	}

	//Constructing the occ table
	//Generate a non-sampled (?) occurence table
	let mut occ_table : HashMap<char, Vec<usize>> = HashMap::new();

	//Populate the occ_table based on the alphabet
	for c in ALPHABET.iter(){
		occ_table.insert(
			*c,
			vec![0; len]
		);
	}

	for i in 0..len{
		let next_occ = &L_col[i];
		//println!("{}", next_occ);

		//Copy the previous row if it isn't the first row
		if i != 0{
			for c in ALPHABET.iter(){
				occ_table.get_mut(c).unwrap()[i] = occ_table[c][i - 1];
			}
		}

		//Update the new row
		//Skip over the terminating character and just leave the copied row
		if *next_occ != '$' {
			occ_table.get_mut(next_occ).unwrap()[i] += 1;
		}
	}

	let mut F_col_compressed = HashMap::<char, usize>::new();
	for c in ALPHABET.iter(){
		F_col_compressed.insert(*c, 0);
	}
	F_col_compressed.insert('$', 0);

	//Compress the first column into a hashmap corresponding to character and index start
	for i in 0..F_col.len(){
		let key = F_col[i];
		*F_col_compressed.get_mut(&key).unwrap() += 1;
	}

	//println!("{:?}", F_col);
	//println!("{:?}", F_col_compressed);

	let fm_index = FmIndex {
		name : header,
		len : len,
		first_column : F_col_compressed,
		last_column : L_col,
		occ : occ_table,
		suffix_arr : suffix_arr_off,
		seq: sequence
	};

	println!("Finished computations, writing to file...");

	if test_seq == None {
		//Write the struct to file
		let encoded: Vec<u8> = bincode::serialize(&fm_index).unwrap();

		let mut file = File::create(output_file).unwrap();
		file.write_all(&encoded).unwrap();
	}

	return fm_index;
}

fn substring(input_file : &str, start : usize, end : usize){
	let mut header : String = String::new();
	let mut sequence = 
			{
				println!("Reading from \'{}\': from {} to {}", input_file, start, end);

				let mut sequence = fs::read_to_string(input_file)
					.expect("Something went wrong reading the file");

				//Read in the sequence and header data from file
				let mut split = sequence.split("\n");

				let mut sequence : String = String::new(); 

				for (i, item) in split.enumerate() {
					if i == 0 { header.push_str(item) }
					else{
						sequence.push_str(item);
					}
				}
				sequence = sequence.trim().to_string();
				//Add the terminating character after processing
				sequence.push_str("$");

				println!("{}", header);
				sequence
		};

	println!("{}", &sequence[start..end]);
}