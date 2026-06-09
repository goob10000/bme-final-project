#r "nuget:Dmx.Amyris.Bio"
#r "nuget: FSharp.Collections.ParallelSeq"

open Amyris.Bio.biolib
open System.Text.RegularExpressions
open System.IO
open FSharp.Collections.ParallelSeq

let root = "c:/proj/Nathaniel/Flu/"
let input = "sequences_partial.fasta"
// let input = "test.fasta"

[<Literal>]
let orfRE = "(ATG((AAA|AAT|AAC|AAG|ATA|ATT|ATC|ATG|ACA|ACT|ACC|ACG|AGA|AGT|AGC|AGG|TAT|TAC|TTA|TTT|TTC|TTG|TCA|TCT|TCC|TCG|TGT|TGC|TGG|CAA|CAT|CAC|CAG|CTA|CTT|CTC|CTG|CCA|CCT|CCC|CCG|CGA|CGT|CGC|CGG|GAA|GAT|GAC|GAG|GTA|GTT|GTC|GTG|GCA|GCT|GCC|GCG|GGA|GGT|GGC|GGG)+)(TAA|TAG|TGA))"

let regexp = Regex(orfRE)
let findORFs (seq: string) =
    let matches = regexp.Matches(seq)
    [ for m in matches -> m.Value ]

let start = System.DateTime.Now
let mutable progress = 0
do 
    use outF = new StreamWriter(Path.Combine(root, "orfs.txt"))
    // process filesand
    fastaStream (root + input) 
    |> PSeq.mapi (fun i (name,seq) ->
        let orfs = findORFs seq
        // if progress % 10000 = 0 then printfn "Processed %d sequences" progress
        // progress <- progress + 1
        
        findORFs seq |> Seq.map (fun orf -> $"{name}\t{orf}") |> String.concat "\n"
    ) 
    |> Seq.iter (fun line -> outF.WriteLine(line))

let elapsed = System.DateTime.Now - start

printfn "Elapsed time: %A" elapsed
printfn $"Processed {progress} sequences in {elapsed.TotalSeconds} seconds ({float progress / elapsed.TotalSeconds} sequences/sec)"