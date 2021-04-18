//! # Assignment 1 Bioinformatica, Davide Cozzi, 829827
//! Crate relativo a parte del primo assignment di Bioinformatica per l'anno 2020/2021.


/// Primo esercizio
pub mod es1 {
    /// Funzione per la creazione dei token.
    /// # Examples
    /// ```
    /// let seq = "aagggcccccttttcc".to_string();
    /// let split = esbio::es1::splitter(&seq);
    /// let check = vec!["aa", "ggg", "ccccc", "tttt","cc"];
    /// assert_eq!(check, split);
    /// ```
    #[allow(dead_code)]
    pub fn splitter(string: &String) -> Vec<String> {
        let chars: Vec<char> = string.chars().collect();
        let mut result = Vec::new();
        let mut tmp = "".to_string();
        if chars.len() == 1 {
            let res = vec![String::from(chars[0])];
            return res;
        }
        for i in 0..chars.len() {
            if i > 0 && chars[i - 1] != chars[i] {
                result.push(tmp.clone());
                tmp = "".to_string();
            }
            tmp.push(chars[i]);
            if i == chars.len() - 1 {
                result.push(tmp.clone());
            }
        }
        result
    }

    /// Funzione per la versione 1 del controllo delle mutazioni, versione basata su token dell'esercizio 1.
    /// # Examples
    /// ```
    /// let seq1 = "ATAGCTC".to_string();
    /// let seq2 = "AAATAAAGGGGCCCCCTTTTTTTCC".to_string();
    /// assert!(esbio::es1::check_virus(&seq1, &seq2));
    /// ```
    /// ```
    /// let seq1 = "TTAGCTC".to_string();
    /// let seq2 = "AAATAAAGGGGCCCCCTTTTTTTCC".to_string();
    /// assert!(!esbio::es1::check_virus(&seq1, &seq2));
    /// ```
    #[allow(dead_code)]
    pub fn check_virus(s1: &String, s2: &String) -> bool {
        let seq1 = s1.to_lowercase();
        let seq2 = s2.to_lowercase();
        let mut check = true;
        if seq1.is_empty() || seq2.is_empty() {
            return false;
        }
        let vecseq1 = splitter(&seq1);
        let vecseq2 = splitter(&seq2);
        if vecseq1.len() != vecseq2.len() {
            return false;
        }
        for (x, y) in vecseq1.iter().zip(vecseq2.iter()) {
            if x.chars().next().unwrap() != y.chars().next().unwrap() || !check {
                return false;
            }
            match x.chars().next().unwrap() {
                'a' => check = check_a(x, y),
                't' => check = check_t(x, y),
                'g' | 'c' => check = check_gc(x, y),
                _ => continue,
            }
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


    /// Funzione per la versione 2 del controllo delle mutazioni, versione basata su indici dell'esercizio 1.
    /// # Examples
    /// ```
    /// let seq1 = "ATAGCTC".to_string();
    /// let seq2 = "AAATAAAGGGGCCCCCTTTTTTTCC".to_string();
    /// assert!(esbio::es1::check_virus2(&seq1, &seq2));
    /// ```
    /// ```
    /// let seq1 = "TTAGCTC".to_string();
    /// let seq2 = "AAATAAAGGGGCCCCCTTTTTTTCC".to_string();
    /// assert!(!esbio::es1::check_virus2(&seq1, &seq2));
    /// ```
    #[allow(dead_code)]
    pub fn check_virus2(s1: &String, s2: &String) -> bool {
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
        for i in 0..seq1.len() {
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
}

/// Secondo esercizio
///
/// Il codice relativo ai grafi di De Bruijn si deve assumere funzionante solo nei casi
/// coperti, tramite le varie assunzioni, dal secondo esercizio
pub mod es2 {
    use std::collections::{HashSet, HashMap};
    use std::fs::File;
    use std::io::Write;

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
        v
    }

    /// Funzione per la creazione del kmer-set a partire da un vettore di stringhe.
    ///
    /// Le sequenze vengono consumate quindi conviene passare cloni.
    /// # Examples
    /// ```
    /// let k = 6;
    /// let seq = "atcttgcat".to_string();
    /// let kmers = vec!["atcttg", "tcttgc", "cttgca", "ttgcat"];
    /// assert_eq!(kmers, esbio::es2::get_kmers_unique(&mut vec![seq.clone()], k));
    /// ```
    /// ```
    /// let seq1 = "atcttga".to_string();
    /// let seq2 = "atcttgc".to_string();
    /// let kmers = vec!["atcttg", "tcttga", "tcttgc"];
    /// assert_eq!(kmers, esbio::es2::get_kmers_unique(&mut vec![seq1.clone(), seq2.clone()], 6));
    /// ```
    #[allow(dead_code)]
    pub fn get_kmers_unique(reads: &mut Vec<String>, k: i32) -> Vec<String> {
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

    /// Struct per rappresentare un grafo di De Bruijn secondo le specifiche dell'esercizio 2.
    ///
    /// Grafo rappresentato tramite:
    /// - insieme dei nodi (come hashmap)
    /// - insieme degli eventuali nodi di partenza (come vettore di etichette)
    /// - insieme degli archi (come vettore di coppie etichetta nodo - etichetta nodo)
    /// - insieme dei nodi adiacenti (come hashmap)
    /// - coppie inizio bubble-fine bubble
    /// - cammino di Eulero
    ///
    /// # Example
    /// ```
    /// let k = 6;
    /// let seq = "atcttgcattaccgccccaatc".to_string();
    /// let dbg = esbio::es2::Dbg::new(vec![seq.clone()], k);
    /// assert_eq!(seq, dbg.get_superstring());
    /// ```
    #[allow(dead_code)]
    pub struct Dbg {
        nodes: HashSet<String>,
        nstart: Vec<String>,
        edges: Vec<(String, String)>,
        node_edge_map: HashMap<String, Vec<String>>,
        bubble: Vec<(String, String)>,
        eulerian_path: Vec<String>,
    }

    #[allow(dead_code)]
    impl Dbg {
        /// Metodo per la creazione di un grafo di De bruijn a partire da un vettore di stringhe
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
            let mut bubble = Vec::new();
            let mut curr = ns[0].clone();
            let mut curr1 = ns[0].clone();
            let mut sb = String::new();
            let mut check_bubble = false;
            loop {
                if !enm.contains_key(&curr) || !enm.contains_key(&curr1) ||
                    enm[&curr].len() == 0 || enm[&curr1].len() == 0 {
                    break;
                } else if enm[&curr].len() == 2 {
                    sb = curr.clone();
                    curr1 = enm[&curr][1].clone();
                    curr = enm[&curr][0].clone();
                    check_bubble = true;
                } else {
                    curr = enm[&curr][0].clone();
                    curr1 = enm[&curr1][0].clone();
                    if check_bubble && curr == curr1 {
                        bubble.push((sb.clone(), curr.clone()));
                        check_bubble = false;
                    }
                }
            }
            let mut dbg = Dbg {
                nodes: n,
                nstart: ns,
                edges: e,
                node_edge_map: enm,
                bubble,
                eulerian_path: vec![],
            };
            let ep = dbg.eulerian_trail();
            dbg.eulerian_path = ep;
            dbg
        }

        /// Metodo per la creazione di un grafo di De bruijn a partire da una singola stringa
        pub fn new_single(seq: &String, k: i32) -> Dbg {
            Dbg::new(vec![seq.to_string()], k)
        }

        /// Metodo per il calcolo del cammino Euleriano
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
        /// Metodo per l'ottenimento della possibile superstringa a partire dal grafo, in presenza di
        /// cammino Euleriano
        pub fn get_superstring(&self) -> String {
            let trail = &self.eulerian_path;
            let mut superstring = trail[0][..trail[0].len() - 1].to_string();
            for elem in trail {
                superstring.push(elem.chars().last().unwrap());
            }
            //println!("{:?}", superstring);
            superstring
        }

        /// Metodo getter per i nodi
        pub fn nodes(&self) -> &HashSet<String> {
            &self.nodes
        }

        /// Metodo getter per i nodi di partenza
        pub fn nstart(&self) -> &Vec<String> {
            &self.nstart
        }

        /// Metodo getter per gli archi
        pub fn edges(&self) -> &Vec<(String, String)> {
            &self.edges
        }

        /// Metodo getter per la mappa dei nodi di adiacebza
        pub fn node_edge_map(&self) -> &HashMap<String, Vec<String>> {
            &self.node_edge_map
        }

        /// Metodo getter per il cammino Euleriano
        pub fn eulerian_path(&self) -> &Vec<String> { &self.eulerian_path }

        /// Metodo getter per ottenere eventyale ordine degli inizi delle bubble
        pub fn bubble(&self) -> &Vec<(String, String)> { &self.bubble }

        /// Metodo per la stampa del grafo in formato DOT
        ///
        /// `output` è l'argomento contentente il path, `bubble` se si vuole la stampa degli archi
        /// che chiudono le bolle
        pub fn to_dot(&self, output: String, bubble: bool) {
            let mut fileout = File::create(output).expect("error");
            fileout
                //.write("digraph sample{\n rankdir=\"LR\";\n".as_bytes())
                .write("digraph sample{\n nodesep=\"0.3\";\nranksep=\"0.3\";\n ".as_bytes())
                .expect("error");
            for edge in &self.edges {
                fileout
                    .write(
                        format!(
                            "\t\"{}\" -> \"{}\" [ label = \" {}\" ];\n",
                            edge.0,
                            edge.1,
                            edge.1.chars().last().unwrap()
                        )
                            .as_bytes(),
                    )
                    .expect("error");
            }
            if bubble {
                for edge in &self.bubble {
                    fileout
                        .write(
                            format!(
                                "\t\"{}\" -> \"{}\" [ style = \" dashed\" ];\n",
                                edge.0,
                                edge.1
                            )
                                .as_bytes(),
                        )
                        .expect("error");
                }
            }
            fileout.write("}".as_bytes()).expect("error");
        }
    }


    /// Struct per rappresentare una mutazione.
    ///
    /// Mutazione rappresentata con:
    /// - base prima sequenza
    /// - base seconda sequenza
    /// - indice della mutazione
    #[derive(Debug)]
    pub struct Mutation {
        ibase: char,
        fbase: char,
        index: usize,
    }


    #[allow(dead_code)]
    impl Mutation {
        /// Metodo di creazione di una mutazione
        pub fn new(ibase: char, fbase: char, index: usize) -> Self {
            Mutation { ibase, fbase, index }
        }

        /// Metodo getter per la prima base della mutazione
        pub fn ibase(&self) -> char { self.ibase }

        /// Metodo getter per la seconda base della mutazione
        pub fn fbase(&self) -> char { self.fbase }

        /// Metodo getter per l'indice della mutazione
        pub fn index(&self) -> usize { self.index }
    }

    /// Funzione per il controllo di mutazione dell'esercizio 2, versione 1.
    ///
    /// Il codice è relativo al confronto diretto tra kmer e assume un k generico, anche se per l'esercizio
    /// è utile avere k=6.
    ///
    /// # Examples
    /// ```
    /// let seq1 = "ATCTTGCATTACCGCCCCAATC".to_string();
    /// let seq2 = "ATCTTACATTACCGTCCCAACC".to_string();
    /// let muts = esbio::es2::check_mutations(&seq1, &seq2, 6);
    /// let check_index = vec![5, 14, 20];
    /// let index: Vec<usize> = muts.0.iter().map(|e| e.index()).collect();
    /// assert_eq!(muts.1, true);
    /// assert_eq!(muts.0.len(), 3);
    /// assert_eq!(check_index, index)
    /// ```
    /// ```
    /// let seq1 = "TATCTTGCATTACCGCCCCAATC".to_string();
    /// let seq2 = "GATCTTACATTACCGTCCCAACC".to_string();
    /// let muts = esbio::es2::check_mutations(&seq1, &seq2, 6);
    /// let check_index = vec![0, 6, 15, 21];
    /// let index: Vec<usize> = muts.0.iter().map(|e| e.index()).collect();
    /// assert_eq!(muts.1, true);
    /// assert_eq!(muts.0.len(), 4);
    /// assert_eq!(check_index, index)
    /// ```
    #[allow(dead_code)]
    pub fn check_mutations(s1: &String, s2: &String, k: i32) -> (Vec<Mutation>, bool) {
        if s1.len() != s2.len() || s1.len() == 0 || s2.len() == 0 {
            return (Vec::new(), false);
        }
        let mut mutations: Vec<Mutation> = Vec::new();
        let seq1 = s1.to_lowercase();
        let seq2 = s2.to_lowercase();
        let kmer1 = get_kmers_unique(&mut vec![seq1.clone()], k);
        let kmer2 = get_kmers_unique(&mut vec![seq2.clone()], k);
        if kmer1.len() != seq1.len() - ((k - 1) as usize) || kmer2.len() != seq2.len() - ((k - 1) as usize) {
            return (Vec::new(), false);
        }
        let mut curr_ind: usize = 0;
        if seq1.chars().next().unwrap() != seq2.chars().next().unwrap() {
            mutations.push(Mutation::new(seq1.chars().next().unwrap(),
                                         seq2.chars().next().unwrap(),
                                         0));
            curr_ind += 1;
        }
        while curr_ind < kmer1.len() {
            if kmer1[curr_ind].chars().last().unwrap() != kmer2[curr_ind].chars().last().unwrap() {
                mutations.push(Mutation::new(kmer1[curr_ind].chars().last().unwrap(),
                                             kmer2[curr_ind].chars().last().unwrap(),
                                             k as usize - 1 + curr_ind));
                curr_ind += k as usize;
            } else {
                curr_ind += 1;
            }
        }
        (mutations, true)
    }

    /// Funzione per il controllo di mutazione dell'esercizio 2, versione 2.
    ///
    /// Il codice è relativo allo studio iterativo del grafo di De Bruijn, studio possibile solo
    /// grazie alle assunzioni fatte. Per l'esercize è utile avere kmer con k=6.
    ///
    /// # Examples
    /// ```
    /// let seq1 = "ATCTTGCATTACCGCCCCAATC".to_string();
    /// let seq2 = "ATCTTACATTACCGTCCCAACC".to_string();
    /// let muts = esbio::es2::check_mutations2(&seq1, &seq2, 6);
    /// let check_index = vec![5, 14, 20];
    /// let index: Vec<usize> = muts.0.iter().map(|e| e.index()).collect();
    /// assert_eq!(muts.1, true);
    /// assert_eq!(muts.0.len(), 3);
    /// assert_eq!(check_index, index)
    /// ```
    /// ```
    /// let seq1 = "TATCTTGCATTACCGCCCCAATC".to_string();
    /// let seq2 = "GATCTTACATTACCGTCCCAACC".to_string();
    /// let muts = esbio::es2::check_mutations2(&seq1, &seq2, 6);
    /// let check_index = vec![0, 6, 15, 21];
    /// let index: Vec<usize> = muts.0.iter().map(|e| e.index()).collect();
    /// assert_eq!(muts.1, true);
    /// assert_eq!(muts.0.len(), 4);
    /// assert_eq!(check_index, index)
    /// ```
    #[allow(dead_code)]
    pub fn check_mutations2(s1: &String, s2: &String, k: i32) -> (Vec<Mutation>, bool) {
        if s1.len() != s2.len() || s1.len() == 0 || s2.len() == 0 {
            return (Vec::new(), false);
        }
        let mut mutations: Vec<Mutation> = Vec::new();
        let seq1 = s1.to_lowercase();
        let seq2 = s2.to_lowercase();

        let kmer1 = get_kmers_unique(&mut vec![seq1.clone()], k);
        let kmer2 = get_kmers_unique(&mut vec![seq2.clone()], k);
        if kmer1.len() != seq1.len() - ((k - 1) as usize) || kmer2.len() != seq2.len() - ((k - 1) as usize) {
            return (Vec::new(), false);
        }

        let dbg = Dbg::new(vec![seq1.clone(), seq2.clone()], k);
        //println!("{:?}", dbg.bubble);
        let mut curr_ind: usize = 0;
        if dbg.nstart().len() == 2 {
            mutations.push(Mutation::new(dbg.nstart()[0].chars().next().unwrap(),
                                         dbg.nstart()[1].chars().next().unwrap(),
                                         0));
        }
        let start = dbg.nstart[0].clone();
        let mut curr = start;
        loop {
            if !dbg.node_edge_map.contains_key(&curr) || dbg.node_edge_map[&curr].len() == 0 {
                break;
            }
            if dbg.node_edge_map[&curr].len() == 2 {
                mutations.push(Mutation::new(dbg.node_edge_map[&curr][0].chars().last().unwrap(),
                                             dbg.node_edge_map[&curr][1].chars().last().unwrap(),
                                             k as usize - 1 + curr_ind));
            }
            curr = dbg.node_edge_map[&curr][0].clone();
            curr_ind += 1;
        }
        (mutations, true)
    }

    /// Funzione per il controllo di mutazione dell'esercizio 2, versione 3.
    ///
    /// Il codice è relativo allo studio iterativo del grafo di De Bruijn con salto delle bubble.
    /// Anche in questo caso lo studio è possibile solo grazie alle assunzioni fatte. Per l'esercizio è utile avere k=6.
    ///
    /// # Examples
    /// ```
    /// let seq1 = "ATCTTGCATTACCGCCCCAATC".to_string();
    /// let seq2 = "ATCTTACATTACCGTCCCAACC".to_string();
    /// let muts = esbio::es2::check_mutations3(&seq1, &seq2, 6);
    /// let check_index = vec![5, 14, 20];
    /// let index: Vec<usize> = muts.0.iter().map(|e| e.index()).collect();
    /// assert_eq!(muts.1, true);
    /// assert_eq!(muts.0.len(), 3);
    /// assert_eq!(check_index, index)
    /// ```
    /// ```
    /// let seq1 = "TATCTTGCATTACCGCCCCAATC".to_string();
    /// let seq2 = "GATCTTACATTACCGTCCCAACC".to_string();
    /// let muts = esbio::es2::check_mutations3(&seq1, &seq2, 6);
    /// let check_index = vec![0, 6, 15, 21];
    /// let index: Vec<usize> = muts.0.iter().map(|e| e.index()).collect();
    /// assert_eq!(muts.1, true);
    /// assert_eq!(muts.0.len(), 4);
    /// assert_eq!(check_index, index)
    /// ```
    #[allow(dead_code)]
    pub fn check_mutations3(s1: &String, s2: &String, k: i32) -> (Vec<Mutation>, bool) {
        if s1.len() != s2.len() || s1.len() == 0 || s2.len() == 0 {
            return (Vec::new(), false);
        }
        let mut mutations: Vec<Mutation> = Vec::new();
        let seq1 = s1.to_lowercase();
        let seq2 = s2.to_lowercase();

        let kmer1 = get_kmers_unique(&mut vec![seq1.clone()], k);
        let kmer2 = get_kmers_unique(&mut vec![seq2.clone()], k);
        if kmer1.len() != seq1.len() - ((k - 1) as usize) || kmer2.len() != seq2.len() - ((k - 1) as usize) {
            return (Vec::new(), false);
        }

        let dbg = Dbg::new(vec![seq1.clone(), seq2.clone()], k);
        let mut curr_ind: usize = 0;
        let mut bubble_ind = 0;
        if dbg.nstart().len() == 2 {
            mutations.push(Mutation::new(dbg.nstart()[0].chars().next().unwrap(),
                                         dbg.nstart()[1].chars().next().unwrap(),
                                         0));
        }
        let start = dbg.nstart[0].clone();
        let mut curr = start;
        loop {
            if !dbg.node_edge_map.contains_key(&curr) || dbg.node_edge_map[&curr].len() == 0 {
                break;
            }
            if dbg.node_edge_map[&curr].len() == 2 {
                mutations.push(Mutation::new(dbg.node_edge_map[&curr][0].chars().last().unwrap(),
                                             dbg.node_edge_map[&curr][1].chars().last().unwrap(),
                                             k as usize - 1 + curr_ind));

                if bubble_ind >= dbg.bubble.len() {
                    break;
                }
                curr = dbg.bubble[bubble_ind].1.clone();
                curr_ind += k as usize;
                bubble_ind += 1;
            } else {
                curr = dbg.node_edge_map[&curr][0].clone();
                curr_ind += 1;
            }
        }
        (mutations, true)
    }
}

/// Terzo esercizio
pub mod es3 {
    use std::cmp;

    /// Funzione per il controllo della distanza di Hamming tramite problema di Manhattan dell'esercizio 3
    ///
    /// In input:
    /// - prima sequenza
    /// - seconda sequenza
    /// - peso arco orizzontale
    /// - peso arco verticale
    /// - peso arco diagonale
    /// # Examples
    /// ```
    /// let seq1 = "AAACTTTTTTTTT".to_string();
    /// let seq2 = "AAAGTTCCGAGCG".to_string();
    /// assert_eq!(8, esbio::es3::check_manhattan(&seq1, &seq2, 1, 0, 0));
    /// assert_eq!(8, esbio::es3::check_manhattan(&seq1, &seq2, 0, 1, 0));
    /// ```
    #[allow(dead_code)]
    pub fn check_manhattan(s1: &String, s2: &String, h: i32, v: i32, d: i32) -> i32 {
        let seq1 = s1.to_lowercase();
        let seq2 = s2.to_lowercase();
        //let mut distance = 0;
        let m = seq1.len();
        let n = seq2.len();
        let mut mm = vec![vec![0; n + 1]; m + 1];
        let mut b = vec![vec!["-"; n + 1]; m + 1];
        for i in 1..(m + 1) {
            mm[i][0] = mm[i - 1][0] + v;
            b[i][0] = "v";
        }
        for j in 1..(n + 1) {
            mm[0][j] = mm[0][j - 1] + h;
            b[0][j] = "h";
        }
        for i in 1..(m + 1) {
            for j in 1..(n + 1) {
                let tmp = cmp::min(mm[i - 1][j] + v, mm[i][j - 1] + h);
                if seq1.chars().nth(i - 1).unwrap() == seq2.chars().nth(j - 1).unwrap() {
                    mm[i][j] = cmp::min(tmp, mm[i - 1][j - 1] + d);
                    b[i][j] = "d";
                } else {
                    mm[i][j] = tmp;
                }
            }
        }
        let distance = mm[m][n];
        /*for m in M {
            println!("{:?}", m);
        }
        for m in B {
            println!("{:?}", m);
        }*/
        distance
    }
}

/// Quarto esercizio
pub mod es4 {
    use std::cmp;

    /// Funzione per il calcolo della longest common substring tra due stringhe dell'esercizio 4
    ///
    /// # Examples
    /// ```
    /// let seq1 = "AAACGCGCTTTTTCCC".to_string();
    /// let seq2 = "AAAGGGGCGCGCTTTTTAAA".to_string();
    /// assert_eq!("CGCGCTTTTT", esbio::es4::lcs(&seq1, &seq2));
    /// ```
    #[allow(dead_code)]
    pub fn lcs(seq1: &String, seq2: &String) -> String {
        let m = seq1.len();
        let n = seq2.len();
        let mut mm = vec![vec![0; n + 1]; m + 1];
        let mut b = vec![vec!["-"; n + 1]; m + 1];
        let d = 1;
        let h = -1;
        let v = -1;
        //println!("{:?}", M);
        let mut max = 0;
        let mut row = 0;
        let mut col = 0;
        for i in 1..(m + 1) {
            for j in 1..(n + 1) {
                let mut tmp = cmp::max(mm[i - 1][j] + (i as i32) * v, mm[i][j - 1] + (j as i32) * h);
                tmp = cmp::max(tmp, 0);
                if seq1.chars().nth(i - 1).unwrap() == seq2.chars().nth(j - 1).unwrap() {
                    mm[i][j] = cmp::max(tmp, mm[i - 1][j - 1] + d);
                    b[i][j] = "d";
                    if max < mm[i][j] {
                        max = mm[i][j];
                        row = i;
                        col = j;
                    }
                } else {
                    /*if b[i-1][j-1] == "d" || b[i][j-1] == "d" || b[i-1][j] == "d" {
                    M[i][j] = 0;
                }else{
                    M[i][j] = cmp::max(tmp,0);
                    //M[i][j] = tmp;
                }*/
                    mm[i][j] = tmp;
                }
            }
        }
        /*println!("{}", mm[m][n]);
        for row in &mm {
            println!("{:?}", row);
        }
        for row in mm {
            for r in row {
                print!("{}  ", r);
            }
            println!("\n");
            //println!("{:?}", row);
        }
        for row in &b {
            println!("{:?}", row);
        }
        for row in &G {
            println!("{:?}", row);
        }*/
        let mut seq = "".to_string();

        let (mut i, mut j) = (row, col);
        let mut check = b[row][col];
        //println!("check: {} -> {}", check, seq1.chars().nth(i - 1).unwrap());
        while b[i][j] != "-" {
            if check == "d" {
                seq.push(seq1.chars().nth(i - 1).unwrap());
                i -= 1;
                j -= 1;
                check = b[i][j];
            } else {
                continue;
            }
        }
        seq = seq.chars().rev().collect();
        /*println!("{}", seq);
        println!("{}", seq.len());
        println!("{} at [{},{}]", max, row, col);*/
        seq
    }
}

#[cfg(test)]
mod tests {
    use crate::es1::{check_virus,
                     check_virus2};

    use crate::es2::{check_mutations,
                     check_mutations2,
                     check_mutations3};

    use crate::es3::check_manhattan;

    use crate::es4::lcs;

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

    #[test]
    fn test_es2_1() {
        let seq1 = "ATCTTGCATTACCGCCCCAATC".to_string();
        let seq2 = "ATCTTACATTACCGTCCCAACC".to_string();
        let muts = check_mutations(&seq1, &seq2, 6);
        println!("mutations: {:?}", muts);
        let check_index = vec![5, 14, 20];
        let index: Vec<usize> = muts.0.iter().map(|e| e.index()).collect();
        assert_eq!(muts.1, true);
        assert_eq!(muts.0.len(), 3);
        assert_eq!(check_index, index);
    }

    #[test]
    fn test_es2_2() {
        let seq1 = "TATCTTGCATTACCGCCCCAATC".to_string();
        let seq2 = "GATCTTACATTACCGTCCCAACC".to_string();
        let muts = check_mutations(&seq1, &seq2, 6);
        println!("mutations: {:?}", muts);
        let check_index = vec![0, 6, 15, 21];
        let index: Vec<usize> = muts.0.iter().map(|e| e.index()).collect();
        assert_eq!(muts.1, true);
        assert_eq!(muts.0.len(), 4);
        assert_eq!(check_index, index);
    }

    #[test]
    fn test_es2_3() {
        let seq1 = "TATCTTGCATTACCGCCCCAAC".to_string();
        let seq2 = "GATCTTACATTACCGTCCCAACC".to_string();
        let muts = check_mutations(&seq1, &seq2, 5);
        assert_eq!(muts.1, false);
        println!("mutations: {:?}", muts);
        assert_eq!(muts.0.len(), 0);
    }

    #[test]
    fn test_es2_4() {
        let seq1 = "TATCTTGCATTACCGCCCCCAAC".to_string();
        let seq2 = "".to_string();
        let muts = check_mutations(&seq1, &seq2, 5);
        assert_eq!(muts.1, false);
        println!("mutations: {:?}", muts);
        assert_eq!(muts.0.len(), 0);
    }

    #[test]
    fn test_es2_5() {
        let seq1 = "TATCTTGCATTACCGCCCCAAC".to_string();
        let seq2 = "TATCTTGCATTACCGCCCCAAC".to_string();
        let muts = check_mutations(&seq1, &seq2, 5);
        assert_eq!(muts.1, true);
        println!("mutations: {:?}", muts);
        assert_eq!(muts.0.len(), 0);
    }

    #[test]
    fn test_es2_6() {
        let seq1 = "TATCTTGCATTACCGCCCCCAAC".to_string();
        let seq2 = "AATCTTGCATTACCGCCCCCAAA".to_string();
        let muts = check_mutations(&seq1, &seq2, 5);
        println!("mutations: {:?}", muts);
        assert_eq!(muts.1, true);
        assert_eq!(muts.0.len(), 2);
    }

    #[test]
    fn test_es21_1() {
        let seq1 = "ATCTTGCATTACCGCCCCAATC".to_string();
        let seq2 = "ATCTTACATTACCGTCCCAACC".to_string();
        let muts = check_mutations2(&seq1, &seq2, 6);
        println!("mutations: {:?}", muts);
        let check_index = vec![5, 14, 20];
        let index: Vec<usize> = muts.0.iter().map(|e| e.index()).collect();
        assert_eq!(muts.1, true);
        assert_eq!(muts.0.len(), 3);
        assert_eq!(check_index, index);
    }

    #[test]
    fn test_es21_2() {
        let seq1 = "TATCTTGCATTACCGCCCCAATC".to_string();
        let seq2 = "GATCTTACATTACCGTCCCAACC".to_string();
        let muts = check_mutations2(&seq1, &seq2, 6);
        println!("mutations: {:?}", muts);
        let check_index = vec![0, 6, 15, 21];
        let index: Vec<usize> = muts.0.iter().map(|e| e.index()).collect();
        assert_eq!(muts.1, true);
        assert_eq!(muts.0.len(), 4);
        assert_eq!(check_index, index);
    }

    #[test]
    fn test_es21_3() {
        let seq1 = "TATCTTGCATTACCGCCCCAAC".to_string();
        let seq2 = "GATCTTACATTACCGTCCCAACC".to_string();
        let muts = check_mutations2(&seq1, &seq2, 5);
        assert_eq!(muts.1, false);
        println!("mutations: {:?}", muts);
        assert_eq!(muts.0.len(), 0);
    }

    #[test]
    fn test_es21_4() {
        let seq1 = "TATCTTGCATTACCGCCCCCAAC".to_string();
        let seq2 = "".to_string();
        let muts = check_mutations2(&seq1, &seq2, 5);
        assert_eq!(muts.1, false);
        println!("mutations: {:?}", muts);
        assert_eq!(muts.0.len(), 0);
    }

    #[test]
    fn test_es21_5() {
        let seq1 = "TATCTTGCATTACCGCCCCAAC".to_string();
        let seq2 = "TATCTTGCATTACCGCCCCAAC".to_string();
        let muts = check_mutations2(&seq1, &seq2, 5);
        assert_eq!(muts.1, true);
        println!("mutations: {:?}", muts);
        assert_eq!(muts.0.len(), 0);
    }

    #[test]
    fn test_es21_6() {
        let seq1 = "TATCTTGCATTACCGCCCCAAC".to_string();
        let seq2 = "AATCTTGCATTACCGCCCCAAA".to_string();
        let muts = check_mutations2(&seq1, &seq2, 5);
        println!("mutations: {:?}", muts);
        assert_eq!(muts.1, true);
        assert_eq!(muts.0.len(), 2);
    }


    #[test]
    fn test_es22_1() {
        let seq1 = "ATCTTGCATTACCGCCCCAATC".to_string();
        let seq2 = "ATCTTACATTACCGTCCCAACC".to_string();
        let muts = check_mutations3(&seq1, &seq2, 6);
        println!("mutations: {:?}", muts);
        let check_index = vec![5, 14, 20];
        let index: Vec<usize> = muts.0.iter().map(|e| e.index()).collect();
        assert_eq!(muts.1, true);
        assert_eq!(muts.0.len(), 3);
        assert_eq!(check_index, index);
    }

    #[test]
    fn test_es22_2() {
        let seq1 = "TATCTTGCATTACCGCCCCAATC".to_string();
        let seq2 = "GATCTTACATTACCGTCCCAACC".to_string();
        let muts = check_mutations3(&seq1, &seq2, 6);
        println!("mutations: {:?}", muts);
        let check_index = vec![0, 6, 15, 21];
        let index: Vec<usize> = muts.0.iter().map(|e| e.index()).collect();
        assert_eq!(muts.1, true);
        assert_eq!(muts.0.len(), 4);
        assert_eq!(check_index, index);
    }

    #[test]
    fn test_es22_3() {
        let seq1 = "TATCTTGCATTACCGCCCCAAC".to_string();
        let seq2 = "GATCTTACATTACCGTCCCAACC".to_string();
        let muts = check_mutations3(&seq1, &seq2, 5);
        assert_eq!(muts.1, false);
        println!("mutations: {:?}", muts);
        assert_eq!(muts.0.len(), 0);
    }

    #[test]
    fn test_es22_4() {
        let seq1 = "TATCTTGCATTACCGCCCCCAAC".to_string();
        let seq2 = "".to_string();
        let muts = check_mutations3(&seq1, &seq2, 5);
        assert_eq!(muts.1, false);
        println!("mutations: {:?}", muts);
        assert_eq!(muts.0.len(), 0);
    }

    #[test]
    fn test_es22_5() {
        let seq1 = "TATCTTGCATTACCGCCCCAAC".to_string();
        let seq2 = "TATCTTGCATTACCGCCCCAAC".to_string();
        let muts = check_mutations3(&seq1, &seq2, 5);
        assert_eq!(muts.1, true);
        println!("mutations: {:?}", muts);
        assert_eq!(muts.0.len(), 0);
    }

    #[test]
    fn test_es22_6() {
        let seq1 = "TATCTTGCATTACCGCCCCAAC".to_string();
        let seq2 = "AATCTTGCATTACCGCCCCAAA".to_string();
        let muts = check_mutations3(&seq1, &seq2, 5);
        println!("mutations: {:?}", muts);
        assert_eq!(muts.1, true);
        assert_eq!(muts.0.len(), 2);
    }

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
    fn test_es3_3() {
        let seq1 = "AAACTT".to_string();
        let seq2 = "AAACTT".to_string();
        let distance = check_manhattan(&seq1, &seq2, 1, 0, 0);
        println!("hamming: {:?}", distance);
        assert_eq!(distance, 0);
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
        assert_eq!(subs, "CGCGCTTTTT");
    }

    #[test]
    fn test_es4_3() {
        let seq1 = "AAACGCGCTTTTTCCCATCGCGCTTTTTTAAA".to_string();
        let seq2 = "AAAGGGGCGCGCTTTTTAAACGCGCTTTTTTC".to_string();
        let subs = lcs(&seq1, &seq2);
        println!("lcs: {}", subs);
        assert_eq!(subs, "AAACGCGCTTTTT");
    }
}
