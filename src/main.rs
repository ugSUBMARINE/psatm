use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;

/// parse pdb file and calculate pseudoatoms (centroids of catalytic atoms) for each residue - will
/// also calculate the centroid if not all catalytic atoms were found
/// will print info if no centroid for a residue could be calculated (if none of the catalytic
/// atoms are present)
/// :parameter
///     *   file_path: path to the pdb file that should be parsed and converted
/// :return
///     *   res_data: lines for the pdb file with the pseudoatoms
fn pdb_parser(file_path: String) -> Vec<String> {
    // map with catalytic atoms
    let cat_map = get_pa_map();
    // new pseudoatom lines
    let mut res_data: Vec<String> = Vec::new();
    // to test whether new residue is encountered
    let mut prev_line: String = "".to_string();
    let mut prev_hash: String = "".to_string();
    // intermediate storage for atom coordinates of current residue in loop
    let mut x_coords: Vec<f32> = Vec::new();
    let mut y_coords: Vec<f32> = Vec::new();
    let mut z_coords: Vec<f32> = Vec::new();

    // read file and iterate over each line
    let file =
        File::open(&file_path).unwrap_or_else(|_| panic!("File [ {} ] not found", &file_path));
    let buf = BufReader::new(file).lines();
    for (cl, line) in buf.flatten().enumerate() {
        if line.starts_with("ATOM  ") {
            // which residue we are at
            let res_hash = &line[17..26];
            // update res_data if we encounter a new residue
            if res_hash != prev_hash {
                // only if we already have data
                if !prev_line.is_empty() {
                    // number of catalytic atoms found for that residue
                    let num_atoms: f32 = x_coords.len() as f32;
                    if num_atoms > 0.0 {
                        // calculate the centroid of the residues catalytic atoms
                        let mut x_mean: f32 = x_coords.iter().sum();
                        let mut y_mean: f32 = y_coords.iter().sum();
                        let mut z_mean: f32 = z_coords.iter().sum();
                        x_mean /= &num_atoms;
                        y_mean /= &num_atoms;
                        z_mean /= &num_atoms;
                        // format the new coordinates to fit the pdb format
                        let coord_string =
                            format!("{:>8.3}{:>8.3}{:>8.3}", &x_mean, &y_mean, &z_mean);
                        // (pseudo)atom serial number
                        let atom_serial = format!("{:>5}", res_data.len() + 1);
                        // replace original data with pseudoatom data
                        let mut new_line = prev_line.to_string();
                        new_line.replace_range(30..54, &coord_string);
                        new_line.replace_range(6..11, &atom_serial);
                        new_line.replace_range(12..16, " X  ");
                        new_line.replace_range(76..78, " X  ");
                        res_data.push(new_line);
                        // prepare coordinate storage for next residue
                        x_coords.clear();
                        y_coords.clear();
                        z_coords.clear();
                    } else {
                        // if no catalytic atoms were found
                        println!("Not able to calculate pseudoatom for [ {} ]", &prev_hash);
                    }
                }
                prev_line = line.to_string();
                prev_hash = res_hash.to_string();
            }
            // test if current atom of current residue is a catalytic atom and if so add it to the
            // coordinate storage
            let cur_res = line[17..20].replace(' ', "").to_string();
            if let Some(catms) = cat_map.get(&cur_res) {
                let cur_atm = line[12..16].replace(' ', "").to_string();
                if catms.contains(&cur_atm) {
                    let x: &f32 = &line[30..38].replace(' ', "").parse().unwrap_or_else(|_| {
                        panic!("Cannot convert coodinate X from line [ {} ] to f32", cl)
                    });
                    let y: &f32 = &line[38..46].replace(' ', "").parse().unwrap_or_else(|_| {
                        panic!("Cannot convert coodinate Y from line [ {} ] to f32", cl)
                    });
                    let z: &f32 = &line[46..54].replace(' ', "").parse().unwrap_or_else(|_| {
                        panic!("Cannot convert coodinate Z from line [ {} ] to f32", cl)
                    });
                    x_coords.push(*x);
                    y_coords.push(*y);
                    z_coords.push(*z);
                }
            };
        }
    }
    if x_coords.len() > 0 {
        let num_atoms = x_coords.len() as f32;
        // calculate the centroid of the residues catalytic atoms
        let mut x_mean: f32 = x_coords.iter().sum();
        let mut y_mean: f32 = y_coords.iter().sum();
        let mut z_mean: f32 = z_coords.iter().sum();
        x_mean /= &num_atoms;
        y_mean /= &num_atoms;
        z_mean /= &num_atoms;
        // format the new coordinates to fit the pdb format
        let coord_string = format!("{:>8.3}{:>8.3}{:>8.3}", &x_mean, &y_mean, &z_mean);
        // (pseudo)atom serial number
        let atom_serial = format!("{:>5}", res_data.len() + 1);
        // replace original data with pseudoatom data
        let mut new_line = prev_line.to_string();
        new_line.replace_range(30..54, &coord_string);
        new_line.replace_range(6..11, &atom_serial);
        new_line.replace_range(12..16, " X  ");
        new_line.replace_range(76..78, " X  ");
        res_data.push(new_line);
        x_coords.clear();
        y_coords.clear();
        z_coords.clear();
    }
    res_data
}

/// specify which atom per residue is a catalytic atom
/// :parameter
///     *   none
/// :return
///     *   cat_at: map that for RESidue 3letter code to vector of catalytic atoms
fn get_pa_map() -> HashMap<String, Vec<String>> {
    let mut cat_at: HashMap<String, Vec<String>> = HashMap::new();
    cat_at.insert("ALA".to_string(), vec!["CB".to_string()]);
    cat_at.insert("CYS".to_string(), vec!["SG".to_string()]);
    cat_at.insert(
        "ASP".to_string(),
        vec!["OD2".to_string(), "OD1".to_string()],
    );
    cat_at.insert(
        "GLU".to_string(),
        vec!["OE1".to_string(), "OE2".to_string()],
    );
    cat_at.insert(
        "PHE".to_string(),
        vec![
            "CG".to_string(),
            "CD1".to_string(),
            "CE1".to_string(),
            "CD2".to_string(),
            "CE2".to_string(),
            "CZ".to_string(),
        ],
    );
    cat_at.insert("GLY".to_string(), vec!["CA".to_string()]);
    cat_at.insert(
        "HIS".to_string(),
        vec![
            "CG".to_string(),
            "ND1".to_string(),
            "CE1".to_string(),
            "NE2".to_string(),
            "CD2".to_string(),
        ],
    );
    cat_at.insert(
        "ILE".to_string(),
        vec![
            "CB".to_string(),
            "CG1".to_string(),
            "CD1".to_string(),
            "CG2".to_string(),
        ],
    );
    cat_at.insert("LYS".to_string(), vec!["NZ".to_string()]);
    cat_at.insert(
        "LEU".to_string(),
        vec![
            "CB".to_string(),
            "CG".to_string(),
            "CD1".to_string(),
            "CD2".to_string(),
        ],
    );
    cat_at.insert("MET".to_string(), vec!["SD".to_string(), "CE".to_string()]);
    cat_at.insert(
        "ASN".to_string(),
        vec!["OD1".to_string(), "ND2".to_string()],
    );
    cat_at.insert(
        "PRO".to_string(),
        vec![
            "N".to_string(),
            "CA".to_string(),
            "CB".to_string(),
            "CG".to_string(),
            "CD".to_string(),
        ],
    );
    cat_at.insert(
        "GLN".to_string(),
        vec!["OE1".to_string(), "NE2".to_string()],
    );
    cat_at.insert(
        "ARG".to_string(),
        vec![
            "NE".to_string(),
            "CZ".to_string(),
            "NH1".to_string(),
            "NH2".to_string(),
        ],
    );
    cat_at.insert("SER".to_string(), vec!["OG".to_string()]);
    cat_at.insert("THR".to_string(), vec!["OG1".to_string()]);
    cat_at.insert(
        "VAL".to_string(),
        vec!["CB".to_string(), "CG1".to_string(), "CG2".to_string()],
    );
    cat_at.insert(
        "TRP".to_string(),
        vec![
            "CE3".to_string(),
            "CZ3".to_string(),
            "CH2".to_string(),
            "CZ2".to_string(),
            "CE2".to_string(),
            "NE1".to_string(),
            "CD1".to_string(),
            "CG".to_string(),
            "CD2".to_string(),
        ],
    );
    cat_at.insert("TYR".to_string(), vec!["OH".to_string()]);
    cat_at
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() == 3 {
        let in_path = &args[1];
        let out_path = &args[2];
        // get data
        let data = pdb_parser(in_path.to_string());
        // write data to file
        let mut f = File::create(out_path)
            .unwrap_or_else(|_| panic!("Unable to create file at [ {} ] ", &out_path));
        for i in data {
            writeln!(f, "{}", i)
                .unwrap_or_else(|_| panic!("Unable to write line [ {} ] to file", &i));
        }
    } else {
        println!(
        "NAME\n\tpsatm - convert residues to pseudoatoms\nSYNOPSIS\n\tpsatm INFILE OUTFILE\nDESCRIPTION\n\tINFILE is the pdbfile containing the data and OUTFILE is the file where the new data should be stored.\n"
)
    }
}
