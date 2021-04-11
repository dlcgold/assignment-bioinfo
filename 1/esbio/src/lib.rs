use std::collections::{HashSet, HashMap};
use std::fs::File;
use std::io::Write;
use std::cmp::{max, min};
use std::cmp;

#[allow(dead_code)]
fn splitter(string: &String) -> Vec<String> {
    let chars: Vec<char> = string.chars().collect();
    let mut result = Vec::new();
    let mut last_mismatch = 0;
    //println!("{}", string);
    let mut tmp = "".to_string();

    if chars.len() == 1 {
        let res = vec![String::from(chars[0])];
        return res;
    }
    for i in 0..chars.len() {
        //println!("{}",i);

        if i > 0 && chars[i - 1] != chars[i] {
            /*let temp_result: String = chars[last_mismatch..i].iter().collect();
            result.push(temp_result);*/
            result.push(tmp.clone());
            tmp = "".to_string();
            //last_mismatch = i;
        }
        tmp.push(chars[i]);
        if i == chars.len() - 1 {
            /*let temp_result: String = chars[last_mismatch..chars.len()].iter().collect();
            result.push(temp_result);*/
            result.push(tmp.clone());
        }
    }
    result
}

#[allow(dead_code)]
fn check_virus(s1: &String, s2: &String) -> bool {
    let seq1 = s1.to_lowercase();
    let seq2 = s2.to_lowercase();
    let mut check = true;
    if seq1.is_empty() || seq2.is_empty() {
        return false;
    }
    let vecseq1 = splitter(&seq1);
    let vecseq2 = splitter(&seq2);
    //println!("{:?}\n{:?}", vecseq1, vecseq2);
    if vecseq1.len() != vecseq2.len() {
        return false;
    }
    for (x, y) in vecseq1.iter().zip(vecseq2.iter()) {
        //println!("compare {} vs {}", x, y);
        if x.chars().next().unwrap() != y.chars().next().unwrap() || !check {
            return false;
        }
        match x.chars().next().unwrap() {
            'a' => check = check_a(x, y),
            't' => check = check_t(x, y),
            'g' | 'c' => check = check_gc(x, y),
            _ => continue,
        }
        //println!("{}", check);
    }
    check
}

#[allow(dead_code)]
fn check_a(s1: &String, s2: &String) -> bool {
    s2.len() <= 5 * s1.len()
}

#[allow(dead_code)]
fn check_t(s1: &String, s2: &String) -> bool {
    s2.len() <= 10 * s1.len()
}

#[allow(dead_code)]
fn check_gc(s1: &String, s2: &String) -> bool {
    s2.len() >= s1.len()
}

#[derive(Debug)]
struct Mutation {
    ibase: char,
    fbase: char,
    index: usize,
}

impl Mutation {
    pub fn new(ibase: char, fbase: char, index: usize) -> Self {
        Mutation { ibase, fbase, index }
    }
}

#[allow(dead_code)]
fn create_kmers(s: String, k: i32) -> Vec<String> {
    get_kmers(&mut vec![s], k)
}

#[allow(dead_code)]
fn get_kmers(reads: &mut Vec<String>, k: i32) -> Vec<String> {
    let mut v: Vec<String> = Vec::new();
    for s in reads {
        while s.len() >= k as usize {
            v.push(s[0..k as usize].to_string());
            *s = String::from(&s[1..]);
        }
    }
    /*let mut set: HashSet<String> = v.into_iter().collect();
    v = set.into_iter().collect();*/
    v
}

#[allow(dead_code)]
fn get_kmers_unique(reads: &mut Vec<String>, k: i32) -> Vec<String> {
    let mut v: Vec<String> = Vec::new();
    for s in reads {
        while s.len() >= k as usize {
            if !v.contains(&s[0..k as usize].to_string()) {
                v.push(s[0..k as usize].to_string());
            }
            *s = String::from(&s[1..]);
        }
    }
    /*let mut set: HashSet<String> = v.into_iter().collect();
    v = set.into_iter().collect();*/
    v
}

#[derive(Debug)]
struct Kmer {
    mer: String,
    string_prov: HashSet<usize>,
    index: Vec<usize>,
}

impl Kmer {
    pub fn new(mer: String, string_prov: HashSet<usize>, index: Vec<usize>) -> Self {
        Kmer { mer, string_prov, index }
    }
}

fn contains_kmer(v: &Vec<Kmer>, kmer: &String) -> (bool, i32) {
    let check = false;
    for (i, elem) in v.iter().enumerate() {
        if &elem.mer == kmer {
            return (true, i as i32);
        }
    }
    (check, -1)
}

fn get_kmers_unique_map(reads: &mut Vec<String>, k: i32) -> HashMap<String, HashMap<usize, Vec<usize>>> {
    //let mut v: Vec<String> = Vec::new();
    //let mut vkmer = Vec::new();
    let mut index = 0;
    let mut index_str = 0;
    let mut kmermap: HashMap<String, HashMap<usize, Vec<usize>>> = HashMap::new();
    for s in reads {
        while s.len() >= k as usize {
            /*if !v.contains(&s[0..k as usize].to_string()) {
                v.push(s[0..k as usize].to_string());
            }*/
            /*let check = contains_kmer(&vkmer, &s[0..k as usize].to_string());
            if !check.0 {
                let mut set = HashSet::new();
                set.insert(index_str);
                vkmer.push(Kmer::new(s[0..k as usize].to_string(), set, vec![index]))
            }else{
                vkmer[check.1 as usize].string_prov.insert(index_str);
                vkmer[check.1 as usize].index.push(index);
            }*/
            if !kmermap.contains_key(&s[0..k as usize].to_string()) {
                let mut tmp = HashMap::new();
                tmp.insert(index_str, vec![index]);
                kmermap.insert(s[0..k as usize].to_string(), tmp);
            } else {
                let tmp = kmermap.get_mut(&s[0..k as usize].to_string()).unwrap();
                if tmp.contains_key(&index_str) {
                    tmp.get_mut(&index_str).unwrap().push(index);
                } else {
                    tmp.insert(index_str, vec![index]);
                }
            }
            *s = String::from(&s[1..]);
            index += 1;
        }
        index_str += 1;
        index = 0;
    }
    /*let mut set: HashSet<String> = v.into_iter().collect();
    v = set.into_iter().collect();*/
    kmermap
}


#[allow(dead_code)]
fn get_kmers_unique_five(reads: &mut Vec<String>, k: i32) -> Vec<String> {
    let mut v: Vec<String> = Vec::new();
    for s in reads {
        while s.len() >= k as usize {
            if !v.contains(&s[0..k as usize].to_string()) {
                v.push(s[0..k as usize].to_string());
            }
            *s = String::from(&s[6..]);
        }
    }
    /*let mut set: HashSet<String> = v.into_iter().collect();
    v = set.into_iter().collect();*/
    v
}

struct Dbg {
    nodes: HashSet<String>,
    nstart: Vec<String>,
    edges: Vec<(String, String)>,
    node_edge_map: HashMap<String, Vec<String>>,
    eulerian_path: Vec<String>,
}

#[allow(dead_code)]
impl Dbg {
    pub fn new(reads: Vec<String>, k: i32) -> Dbg {
        let kmers = get_kmers_unique(&mut reads.clone(), k);
        let mut n: HashSet<String> = HashSet::new();
        let mut nst: HashSet<String> = HashSet::new();
        let mut e: Vec<(String, String)> = Vec::new();
        for kmer in kmers {
            let kmer1 = &kmer[..kmer.len() - 1];
            let kmer2 = &kmer[1..];
            n.insert(kmer1.to_string());
            n.insert(kmer2.to_string());
            //if !e.contains(&(kmer1.to_string(), kmer2.to_string())) {
            e.push((kmer1.to_string(), kmer2.to_string()));
            //}
            nst.insert(kmer2.to_string());
        }
        let mut ns = Vec::new();
        for node in &n {
            if !nst.contains(node) {
                ns.push(node.to_string());
            }
        }
        let mut enm: HashMap<String, Vec<String>> = HashMap::new();
        for (start, end) in &e {
            if !enm.contains_key(start) {
                enm.insert(start.to_string(), vec![end.to_string()]);
            } else {
                enm.entry(start.to_string())
                    .or_insert(vec![])
                    .push(end.to_string());
            }
        }

        let mut dbg = Dbg {
            nodes: n,
            nstart: ns,
            edges: e,
            node_edge_map: enm,
            eulerian_path: vec![],
        };
        let ep = dbg.eulerian_trail();
        dbg.eulerian_path = ep;
        dbg
    }
    pub fn new_single(seq: &String, k: i32) -> Dbg {
        Dbg::new(vec![seq.to_string()], k)
    }
    pub fn eulerian_trail(&self) -> Vec<String> {
        let mut map = self.node_edge_map.clone();
        let mut trail = Vec::new();
        let mut start = match self.nstart.is_empty() {
            true => self.nodes.iter().next().unwrap().to_string(),
            false => self.nstart.iter().next().unwrap().to_string(),
        };
        trail.push(start.clone());
        loop {
            let mut tmp: Vec<String> = Vec::new();
            let mut prev = start.clone();
            loop {
                if !map.contains_key(&prev) {
                    break;
                }
                let next = map.entry(prev.to_string()).or_insert(vec![]).pop().unwrap();
                if map[&prev.to_string()].is_empty() {
                    map.remove(&prev);
                }
                tmp.push(next.clone());
                if next == start {
                    break;
                }
                prev = next;
            }
            //println!("{:?}", tmp);
            let index = (&trail).iter().position(|r| r == &start).unwrap();
            let mut trailtmp = Vec::new();
            for i in 0..index + 1 {
                trailtmp.push(trail[i].clone());
            }
            for elem in tmp {
                trailtmp.push(elem.clone());
            }
            for i in index + 1..trail.len() {
                trailtmp.push(trail[i].clone());
            }
            trail = trailtmp;
            if map.is_empty() {
                break;
            }
            let mut check_start = false;
            for elem in &trail {
                if map.contains_key(elem) {
                    start = elem.clone();
                    check_start = true;
                    break;
                }
            }
            if !check_start {
                break;
            }
        }
        trail
    }
    pub fn get_superstring(&self) -> String {
        let trail = &self.eulerian_path;
        let mut superstring = trail[0][..trail[0].len() - 1].to_string();
        for elem in trail {
            superstring.push(elem.chars().last().unwrap());
        }
        //println!("{:?}", superstring);
        superstring
    }

    pub fn nodes(&self) -> &HashSet<String> {
        &self.nodes
    }
    pub fn nstart(&self) -> &Vec<String> {
        &self.nstart
    }
    pub fn edges(&self) -> &Vec<(String, String)> {
        &self.edges
    }
    pub fn node_edge_map(&self) -> &HashMap<String, Vec<String>> {
        &self.node_edge_map
    }
    pub fn eulerian_path(&self) -> &Vec<String> {
        &self.eulerian_path
    }
    pub fn to_dot(&self, output: String) {
        let mut fileout = File::create(output).expect("error");
        fileout
            .write("digraph sample{\n".as_bytes())
            .expect("error");
        for edge in &self.edges {
            if edge.0.len() <= 6 {
                fileout
                    .write(
                        format!(
                            "\t\"{}\" -> \"{}\" [ label = \"{}\" ];\n",
                            edge.0,
                            edge.1,
                            edge.1.chars().last().unwrap()
                        )
                            .as_bytes(),
                    )
                    .expect("error");
            } else {
                let mut start = edge.0[0..3].to_string();
                start.push_str(&"...".to_string());
                start.push_str(&edge.0[edge.0.len() - 3..]);
                let mut end = edge.1[0..3].to_string();
                end.push_str(&"...".to_string());
                end.push_str(&edge.1[edge.1.len() - 3..]);
                fileout
                    .write(
                        format!(
                            "\t\"{}\" -> \"{}\" [ label = \"{}\" ];\n",
                            start,
                            end,
                            edge.1.chars().last().unwrap()
                        )
                            .as_bytes(),
                    )
                    .expect("error");
            }
        }
        fileout.write("}".as_bytes()).expect("error");
    }
    pub fn count_bubble(&self) -> u32 {
        let count = 0;
        for e in &self.node_edge_map {
            println!("{:?}", e);
        }
        count
    }
}

fn check_mutations(s1: &String, s2: &String) -> Vec<Mutation> {
    let mut mutations = Vec::new();
    let seq1 = s1.to_lowercase();
    let seq2 = s2.to_lowercase();
    let dbg = Dbg::new(vec![seq1.clone(), seq2.clone()], 5);
    println!("{},{}", &seq1, &seq2);
    //println!("{:?}", dbg.edges);
    //dbg.count_bubble();
    dbg.to_dot("outputs/mut.dot".to_string());

    let mut kmers = get_kmers_unique_map(&mut vec![seq1.clone(), seq2.clone()], 5);
    /*let mut kmers1 = get_kmers_unique(&mut vec![seq1.clone()], 5);
    let mut kmers2 = get_kmers_unique(&mut vec![seq2.clone()], 5);*/
    /*kmers1.sort();
    kmers2.sort();
    kmers.sort();*/
    //println!("kmers1: {:?}\nkmers2: {:?}\nkmers: {:?}", kmers1, kmers2, kmers);
    //println!("kmers1: {:?}\nkmers2: {:?}", kmers1, kmers2);
    /*for kmer in &kmers {
        println!("{:?}", kmer);
    }*/
    //println!("\n");
    /*for e in &dbg.node_edge_map {
        //println!("{:?}", e);

        if (&e.1).len() == 2 {
            //println!("{:?}, {:?}", e.1[0], e.1[1]);
            let vec1 = dbg.node_edge_map.get(&e.1[0]);
            let vec2 = dbg.node_edge_map.get(&e.1[1]);
            let mut tmp1 = e.0.clone();
            tmp1.push(e.1[0].chars().last().unwrap());
            let mut tmp2 = e.0.clone();
            tmp2.push(e.1[1].chars().last().unwrap());
            //println!("{}, {}", &tmp1, &tmp2);
            let index1vec = kmers[&tmp1].get(&0);
            let index2vec = kmers[&tmp2].get(&1);
            //let vec1 = &dbg.node_edge_map[&e.1[0]];
            //println!("{:?},{:?}", index1vec, index2vec);
            let mut good_index = Vec::new();
            if !index1vec.is_none() && !index2vec.is_none() {
                let index1vec = index1vec.unwrap();
                let index2vec = index2vec.unwrap();
                //println!("{:?},{:?}", index1vec, index2vec);
                for ind in index1vec {
                    //if !index2vec.contains(ind) {
                    good_index.push(*ind);
                    //}
                }
                //let mut good_index1 = Vec::new();
                for ind in index2vec {
                    if !index1vec.contains(ind) {
                        good_index.push(*ind);
                    }
                }
                //println!("{:?}", good_index);
                good_index.dedup();
            } else if !index1vec.is_none() && index2vec.is_none() {
                good_index = index1vec.unwrap().to_vec();
            } else {
                good_index = index2vec.unwrap().to_vec();
            }
            println!("{:?}", good_index);
            //let vec2 = &dbg.node_edge_map[&e.1[1]];
            //println!("{:?}", vec2);
            if good_index.contains(&0) {
                let tmp1 = &seq1[0..5].to_string();
                let tmp2 = &seq2[0..5].to_string();
                for (i, x) in tmp1.chars().zip(tmp2.chars()).enumerate() {
                    // println!("{} vs {}", x.1, x.0);
                    if x.0 != x.1 {
                        // println!("insert {},{},{}",x.0, x.1, i);
                        mutations.push(Mutation::new(x.0, x.1, i))
                    }
                }
            } else {
                if vec1.is_none() {
                    //mutations.push(Mutation::new(e.1[0].chars().last().unwrap(), '-', good_index[0] + 4));
                    //println!("insert {},{},{}",e.1[0].chars().last().unwrap(), '-', good_index[0] + 4);
                    if good_index.len() == 2 {
                        mutations.push(Mutation::new('-', e.1[1].chars().last().unwrap(), good_index[1] + 4));
                        //println!("insert {},{},{}",'-', e.1[1].chars().last().unwrap(), good_index[1] + 4);
                    } else {
                        mutations.push(Mutation::new('-', e.1[1].chars().last().unwrap(), good_index[0] + 4));
                        //println!("insert {},{},{}",'-', e.1[1].chars().last().unwrap(), good_index[0] + 4);
                    }
                } else if vec2.is_none() {
                    if good_index.len() == 2 {
                        mutations.push(Mutation::new('-', e.1[1].chars().last().unwrap(), good_index[1] + 4));
                        //println!("insert {},{},{}",'-', e.1[1].chars().last().unwrap(), good_index[1] + 4);
                    } else {
                        mutations.push(Mutation::new('-', e.1[1].chars().last().unwrap(), good_index[0] + 4));
                        //println!("insert {},{},{}",'-', e.1[1].chars().last().unwrap(), good_index[0] + 4);
                    }
                } else {
                    mutations.push(Mutation::new(e.1[0].chars().last().unwrap(), e.1[1].chars().last().unwrap(), good_index[0] + 4));
                    //println!("insert {},{},{}",e.1[0].chars().last().unwrap(), e.1[1].chars().last().unwrap(), good_index[0] + 4);
                }
            }
        }
    }
    /*println!("{:?}", dbg.eulerian_path);
    println!("strt{:?}", dbg.nstart);*/
    if !dbg.nstart.is_empty() {
        for (i, x) in dbg.nstart[0].chars().zip(dbg.nstart[1].chars()).enumerate() {
            if x.0 != x.1 {
                mutations.push(Mutation::new(x.0, x.1, i))
            }
        }
    }*/
    let dbg1 = Dbg::new(vec![seq1.clone()], 4);
    let dbg2 = Dbg::new(vec![seq2.clone()], 4);
    dbg1.to_dot("outputs/mut1.dot".to_string());
    dbg2.to_dot("outputs/mut2.dot".to_string());
    println!("{}\n{}", dbg1.get_superstring(), dbg2.get_superstring());
    mutations
}

fn check_manhattan(s1: &String, s2: &String, h: i32, v: i32, d: i32) -> i32 {
    let seq1 = s1.to_lowercase();
    let seq2 = s2.to_lowercase();
    let mut distance = 0;
    let m = seq1.len();
    let n = seq2.len();
    let mut M = vec![vec![0; n + 1]; m + 1];
    let mut B = vec![vec!["-"; n + 1]; m + 1];
    for i in 1..(m + 1) {
        M[i][0] = M[i - 1][0] + v;
        B[i][0] = "v";
    }
    for j in 1..(n + 1) {
        M[0][j] = M[0][j - 1] + h;
        B[0][j] = "h";
    }
    for i in 1..(m + 1) {
        for j in 1..(n + 1) {
            let tmp = min(M[i - 1][j] + v, M[i][j - 1] + h);
            if seq1.chars().nth(i - 1).unwrap() == seq2.chars().nth(j - 1).unwrap() {
                M[i][j] = min(tmp, M[i - 1][j - 1] + d);
                B[i][j] = "d";
            } else {
                M[i][j] = tmp;
            }
        }
    }
    distance = M[m][n];
    for m in M {
        println!("{:?}", m);
    }
    for m in B {
        println!("{:?}", m);
    }
    distance
}

fn lcs_old(seq1: &String, seq2: &String) -> String {
    /*let seq1 = s1.to_lowercase();
    let seq2 = s2.to_lowercase();*/
    let m = seq1.len();
    let n = seq2.len();
    //println!("{}, {}", n, m);
    let mut M = vec![vec![0; n + 1]; m + 1];
    let mut B = vec![vec!["-"; n + 1]; m + 1];
    let mut G = vec![vec![1; n + 1]; m + 1];

    //println!("{:?}", M);
    let mut max = 0;
    let mut row = 0;
    let mut col = 0;
    for i in 0..(m + 1) {
        for j in 0..(n + 1) {
            if i == 0 || j == 0 {
                M[i][j] = 0;
            } else if seq1.chars().nth(i - 1).unwrap() == seq2.chars().nth(j - 1).unwrap() {
                M[i][j] = M[i - 1][j - 1] + 1;
                B[i][j] = "d";
                G[i][j] = 0;
                if max < M[i][j] {
                    max = M[i][j];
                    row = i;
                    col = j;
                }
            } else {
                M[i][j] = -1;
                /*M[i][j] = max(M[i - 1][j], M[i][j - 1]);
                if M[i][j] == M[i - 1][j] {
                    B[i][j] = "u";
                } else if M[i][j] == M[i][j - 1] {
                    B[i][j] = "l";
                }*/
            }
        }
    }
    println!("{}", M[m][n]);
    for row in M {
        println!("{:?}", row);
    }
    for row in &B {
        println!("{:?}", row);
    }
    /*for row in &G {
        println!("{:?}", row);
    }*/
    let mut seq = "".to_string();

    let (mut i, mut j) = (row, col);
    let mut check = B[row][col];
    println!("check: {} -> {}", check, seq1.chars().nth(i - 1).unwrap());
    while B[i][j] != "-" {
        if check == "d" {
            seq.push(seq1.chars().nth(i - 1).unwrap());
            i -= 1;
            j -= 1;
            check = B[i][j];
        } else {
            continue;
        }
    }
    seq = seq.chars().rev().collect();
    println!("{}", seq);
    println!("{}", seq.len());
    println!("{} at [{},{}]", max, row, col);
    seq
}

#[allow(dead_code)]
fn lcs(seq1: &String, seq2: &String) -> String {
    /*let seq1 = s1.to_lowercase();
    let seq2 = s2.to_lowercase();*/
    let m = seq1.len();
    let n = seq2.len();
    //println!("{}, {}", n, m);
    let mut M = vec![vec![0; n + 1]; m + 1];
    let mut B = vec![vec!["-"; n + 1]; m + 1];
    let mut G = vec![vec![1; n + 1]; m + 1];
    let d = 1;
    let h = -1;
    let v = -1;
    //println!("{:?}", M);
    let mut max = 0;
    let mut row = 0;
    let mut col = 0;
    for i in 1..(m + 1) {
        for j in 1..(n + 1) {
            let mut tmp = cmp::max(M[i - 1][j] + (i as i32) * v, M[i][j - 1] + (j as i32) * h);
            tmp = cmp::max(tmp, 0);
            if seq1.chars().nth(i - 1).unwrap() == seq2.chars().nth(j - 1).unwrap() {
                M[i][j] = cmp::max(tmp, M[i - 1][j - 1] + d);
                B[i][j] = "d";
                if max < M[i][j] {
                    max = M[i][j];
                    row = i;
                    col = j;
                }
            } else {
                /*if B[i-1][j-1] == "d" || B[i][j-1] == "d" || B[i-1][j] == "d" {
                    M[i][j] = 0;
                }else{
                    M[i][j] = cmp::max(tmp,0);
                    //M[i][j] = tmp;
                }*/
                M[i][j] = tmp;
            }
        }
    }
    println!("{}", M[m][n]);
    for row in &M {
        println!("{:?}", row);
    }
    for row in M {
        for r in row {
            print!("{}  ", r);
        }
        println!("\n");
        //println!("{:?}", row);
    }
    for row in &B {
        println!("{:?}", row);
    }
    /*for row in &G {
        println!("{:?}", row);
    }*/
    let mut seq = "".to_string();

    let (mut i, mut j) = (row, col);
    let mut check = B[row][col];
    println!("check: {} -> {}", check, seq1.chars().nth(i - 1).unwrap());
    while B[i][j] != "-" {
        if check == "d" {
            seq.push(seq1.chars().nth(i - 1).unwrap());
            i -= 1;
            j -= 1;
            check = B[i][j];
        } else {
            continue;
        }
    }
    seq = seq.chars().rev().collect();
    println!("{}", seq);
    println!("{}", seq.len());
    println!("{} at [{},{}]", max, row, col);
    seq
}

#[allow(dead_code)]
fn check_virus2(s1: &String, s2: &String) -> bool {
    let seq1 = s1.to_lowercase();
    let seq2 = s2.to_lowercase();
    let mut check = true;
    if seq1.is_empty() || seq2.is_empty() {
        return false;
    }
    let mut count1 = 0;
    let mut count2 = 0;
    let seq1: Vec<char> = seq1.chars().collect();
    let seq2: Vec<char> = seq2.chars().collect();
    let mut j = 0;
    for mut i in 0..seq1.len() {
        //println!("{}", i);
        count1 += 1;
        //println!("s1: {}{}", seq1[i], seq1[i + 1]);
        /* if i == seq1.len() - 1 {
             println!("{},{}|{},{}", i, seq1.len(), j, seq2.len());
             if i == seq1.len() && j != seq2.len() {
                 return false;
             }
         }*/
        if i == seq1.len() - 1 || seq1[i] != seq1[i + 1] {
            //println!("s12: {}{}", seq1[i], seq2[j]);
            if seq1[i] != seq2[j] {
                return false;
            }

            //println!("s2: {}{}", seq2[j], seq2[j + 1]);
            //j += 1;
            while j != seq2.len() - 1 && seq2[j] == seq2[j + 1] {
                //println!("s2: {}{}", seq2[j], seq2[j + 1]);
                count2 += 1;
                j += 1;
                /*if j == seq2.len() - 1 {
                    break;
                }*/
            }
            count2 += 1;
            j += 1;

            //println!("{}, {} for {}", count1, count2, seq1[i]);

            if seq1[i] == 'a' && count2 <= 5 * count1 {
                check = true;
            } else if seq1[i] == 't' && count2 <= 10 * count1 {
                check = true
            } else if (seq1[i] == 'c' || seq1[i] == 'g') && count2 >= count1 {
                check = true
            } else {
                return false;
            }
            count2 = 0;
            count1 = 0;
        }
        //println!("{},{}|{},{}", i, seq1.len(), j, seq2.len());
        if i == seq1.len() - 1 && j != seq2.len() {
            return false;
        }
    }
    // println!("check {}", check);
    check
}

#[cfg(test)]
mod tests {
    use crate::{check_virus, check_mutations, get_kmers_unique_map, check_manhattan, lcs, check_virus2, lcs_old};

    #[test]
    fn test_es1_1() {
        let seq1 = "ATAGCTC".to_string();
        let seq2 = "AAATAAAGGGGCCCCCTTTTTTTCC".to_string();
        assert!(check_virus(&seq1, &seq2));
    }

    #[test]
    fn test_es1_2() {
        let seq1 = "TTAGCTC".to_string();
        let seq2 = "AAATAAAGGGGCCCCCTTTTTTTCC".to_string();
        assert!(!check_virus(&seq1, &seq2));
    }

    #[test]
    fn test_es1_3() {
        let seq1 = "ATAGCTC".to_string();
        let seq2 = "AAATAAAGGGGCCCCCTTTTTTTTTTTTTTTCC".to_string();
        assert!(!check_virus(&seq1, &seq2));
    }

    #[test]
    fn test_es1_4() {
        let seq1 = "ATAGCTC".to_string();
        let seq2 = "AAAAAATAAAGGGGCCCCCTTTTTTTCC".to_string();
        assert!(!check_virus(&seq1, &seq2));
    }

    #[test]
    fn test_es1_5() {
        let seq1 = "ATAGCC".to_string();
        let seq2 = "AAATAAAGGGGCCCCCTTTTTTTCC".to_string();
        assert!(!check_virus(&seq1, &seq2));
    }

    #[test]
    fn test_es1_6() {
        let seq1 = "ATAAGCTC".to_string();
        let seq2 = "AAATAAAAAAAAAAGGGGCCCCCTTTTTTTCC".to_string();
        assert!(check_virus(&seq1, &seq2));
    }

    #[test]
    fn test_es1_7() {
        let seq1 = "ATAAGCTC".to_string();
        let seq2 = "AAACCCTAAAAAAAAAAGGGGCCCCCTTTTTTT".to_string();
        assert!(!check_virus(&seq1, &seq2));
    }

    #[test]
    fn test_es1_8() {
        let seq1 = "".to_string();
        let seq2 = "AAATAAAAAAAAAAGGGGCCCCCTTTTTTTCC".to_string();
        assert!(!check_virus(&seq1, &seq2));
    }

    #[test]
    fn test_es1_9() {
        let seq1 = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        let seq2 = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        assert!(check_virus(&seq1, &seq2));
    }

    #[test]
    fn test_es2_1() {
        /*let seq1 = "AATTAAATTAAAAAGCGGCCCGCTTTATTTCC".to_string();
        let seq2 = "AAATAAAAAAAAAAGGGGCCCCCTTTTTTTCC".to_string();*/
        let seq1 = "AATTAAAGT".to_string();
        let seq2 = "ATTTAAAAT".to_string();
        let muts = check_mutations(&seq1, &seq2);
        println!("mutations: {:?}", muts);
        assert_eq!(2, muts.len());
    }

    #[test]
    fn test_es2_2() {
        let seq1 = "AATTAAAGTTTAACC".to_string();
        let seq2 = "ATTTAAAATTTAAC".to_string();
        let muts = check_mutations(&seq1, &seq2);
        println!("mutations: {:?}", muts);
        assert_eq!(3, muts.len());
    }

    /*#[test]
    fn test_es2_3() {
        /*let seq1 = "AATTAAATTAAAAAGCGGCCCGCTTTATTTCC".to_string();
        let seq2 = "AAATAAAAAAAAAAGGGGCCCCCTTTTTTTCC".to_string();*/
       /* let seq1 = "AATTAAAGT".to_string();
        let seq2 = "ATTTAAAAT".to_string();*/
        let seq1 = "AATTAAAGTTTAA".to_string();
        let seq2 = "ATTTAAAATTTAAC".to_string();
        let muts = get_kmers_unique_map(&mut vec![seq1.clone(), seq2.clone()], 5);
        println!("mutations: {:?}", muts);
        for elem in muts{
            print!("{:?}\n", elem);
        }
        let muts = check_mutations(&seq1, &seq2);
        println!("mutations: {:?}", muts);
        assert_eq!(2, 2);
    }*/

    #[test]
    fn test_es3_1() {
        let seq1 = "AAACTTT".to_string();
        let seq2 = "AAAGTTC".to_string();
        let distance = check_manhattan(&seq1, &seq2, 1, 0, 0);
        println!("hamming: {:?}", distance);
        assert_eq!(distance, 2);
    }

    #[test]
    fn test_es3_2() {
        let seq1 = "AAACTTTTTTTTT".to_string();
        let seq2 = "AAAGTTCCGAGCG".to_string();
        let distance = check_manhattan(&seq1, &seq2, 1, 0, 0);
        println!("hamming: {:?}", distance);
        assert_eq!(distance, 8);
    }

    #[test]
    fn test_es4_1() {
        let seq1 = "AAACGCGCTTTTTCCC".to_string();
        let seq2 = "AAAGGGGCGCGCTTTTTAAA".to_string();
        let subs = lcs(&seq1, &seq2);
        println!("lcs: {}", subs);
        assert_eq!(subs, "CGCGCTTTTT");
    }

    #[test]
    fn test_es4_2() {
        let seq1 = "AAACGCGCTTTTTCCCAT".to_string();
        let seq2 = "AAAGGGGCGCGCTTTTTAAA".to_string();
        let subs = lcs(&seq1, &seq2);
        println!("lcs: {}", subs);
        //assert_eq!(3, 3);
        assert_eq!(subs, "CGCGCTTTTT");
    }

    #[test]
    fn test_es4_3() {
        let seq1 = "AAACGCGCTTTTTCCCATCGCGCTTTTTTAAA".to_string();
        let seq2 = "AAAGGGGCGCGCTTTTTAAACGCGCTTTTTTC".to_string();
        let subs = lcs(&seq1, &seq2);
        println!("lcs: {}", subs);
        //assert_eq!(3, 3);
        assert_eq!(subs, "AAACGCGCTTTTT");
    }



    #[test]
    fn test_es11_1() {
        let seq1 = "ATAGCTC".to_string();
        let seq2 = "AAATAAAGGGGCCCCCTTTTTTTCC".to_string();
        assert!(check_virus2(&seq1, &seq2));
    }

    #[test]
    fn test_es11_2() {
        let seq1 = "TTAGCTC".to_string();
        let seq2 = "AAATAAAGGGGCCCCCTTTTTTTCC".to_string();
        assert!(!check_virus2(&seq1, &seq2));
    }

    #[test]
    fn test_es11_3() {
        let seq1 = "ATAGCTC".to_string();
        let seq2 = "AAATAAAGGGGCCCCCTTTTTTTTTTTTTTTCC".to_string();
        assert!(!check_virus2(&seq1, &seq2));
    }

    #[test]
    fn test_es11_4() {
        let seq1 = "ATAGCTC".to_string();
        let seq2 = "AAAAAATAAAGGGGCCCCCTTTTTTTCC".to_string();
        assert!(!check_virus2(&seq1, &seq2));
    }

    #[test]
    fn test_es11_5() {
        let seq1 = "ATAGCC".to_string();
        let seq2 = "AAATAAAGGGGCCCCCTTTTTTTCC".to_string();
        assert!(!check_virus2(&seq1, &seq2));
    }

    #[test]
    fn test_es11_6() {
        let seq1 = "ATAAGCTC".to_string();
        let seq2 = "AAATAAAAAAAAAAGGGGCCCCCTTTTTTTCC".to_string();
        assert!(check_virus2(&seq1, &seq2));
    }

    #[test]
    fn test_es11_7() {
        let seq1 = "ATAAGCTC".to_string();
        let seq2 = "AAACCCTAAAAAAAAAAGGGGCCCCCTTTTTTT".to_string();
        assert!(!check_virus2(&seq1, &seq2));
    }

    #[test]
    fn test_es11_8() {
        let seq1 = "".to_string();
        let seq2 = "AAATAAAAAAAAAAGGGGCCCCCTTTTTTTCC".to_string();
        assert!(!check_virus2(&seq1, &seq2));
    }

    #[test]
    fn test_es11_9() {
        let seq1 = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        let seq2 = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC".to_string();
        assert!(check_virus2(&seq1, &seq2));
    }
}
