import json
import polars as pl
data = json.load(open("Genomes/Covid/ncov_open_global_6m.json"))

data.keys()
data["tree"]["children"][0]["node_attrs"]["epiweek"]
data["tree"]["children"][0]["branch_attrs"]["mutations"]
data["tree"]["children"][1]["children"][1].keys()
data["tree"]["children"][2].keys()

for x in data["tree"]["children"][1]["children"]:
    print(x.keys())
    print(x["name"])


def recurse_wrapper(node:dict):
    node_name = []
    node_mutations = []
    node_epiweek = []
    def merge_mutations(parent_mutations:dict, current_mutations:dict):
        merged = {key: values[:] for key, values in parent_mutations.items()}
        for gene, mutations in current_mutations.items():
            merged.setdefault(gene, [])
            merged[gene].extend(mutations)
        return merged

    def recurse(node_list:list, inherited_mutations:dict):
        # print("recursing")
        for i in node_list:
            # print(f"{i=}")
            # print(f"{i.keys()=}")
            current_mutations = merge_mutations(inherited_mutations, i["branch_attrs"]["mutations"])
            if "children" not in i.keys():
                # print(i.keys())
                node_name.append(i["name"])
                node_mutations.append(current_mutations)
                node_epiweek.append(i["node_attrs"]["epiweek"])
            else:
                recurse(i["children"], current_mutations)
    recurse(node["tree"]["children"], {})
    return node_name, node_mutations, node_epiweek

nn, nm, ne = recurse_wrapper(data)
nm
weeks_since_W1 = [(int(d["value"][0:4]) - 2020)*52 + int(d["value"][4:]) for d in ne]




df = pl.DataFrame({"name":nn, "mutations":nm, "weeks":weeks_since_W1})

df.filter(pl.col("mutations") != df["mutations"][1])
df.write_parquet("Genomes/Covid/mutations2.parquet")