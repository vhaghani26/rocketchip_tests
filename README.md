# Rocketchip Tests

(clone the environment and repository information)

## Expected vs. Observed Peaks

First, set up the sequencing data. To do so, enter `exp_vs_obs/seq_data/`. Run the following to unzip `seqdata.tar.gz`:

```
tar -xvzf seqdata.tar.gz
```

Now run the following to reorganize the data for use:

```
mkdir paired_broad paired_narrow single_broad single_narrow

for dir in paired_broad_*; do
    new_name="test_$(echo "$dir" | sed 's/paired_broad_//')"
    mv "$dir" "paired_broad/$new_name"
done

for dir in paired_narrow_*; do
    new_name="test_$(echo "$dir" | sed 's/paired_narrow_//')"
    mv "$dir" "paired_narrow/$new_name"
done

for dir in single_broad_*; do
    new_name="test_$(echo "$dir" | sed 's/single_broad_//')"
    mv "$dir" "single_broad/$new_name"
done

for dir in single_narrow_*; do
    new_name="test_$(echo "$dir" | sed 's/single_narrow_//')"
    mv "$dir" "single_narrow/$new_name"
done
```

Navigate back to `exp_vs_obs/` and open `exp_vs_obs.py` in a text editor of your choice. In the "Import Modules" section, go to the "User-specific variables" section and change the working directory and authors variable to your working directory and name. Save and close the file.

Now, run the script. This will trigger the creation of all project files, creation of all snakefiles, execute the snakefiles to conduct the analysis, carry out statistics, and visualize the results. Effectively, this will run the test beginning to end. Please not this will take a long time due to the amount of tests:

```
python3 exp_vs_obs.py
```
