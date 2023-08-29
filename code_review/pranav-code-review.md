# Installation

On Pop!\_OS 22.04, which is the same as Ubuntu 22.04:
    - First ran `sudo apt install libgdal-dev`, without which `rgdal` would not
      install.
    - Remaining `install.packages(...)` as per the README provided, worked
      perfectly.

# Functions

Below, I have gone through the functions file-by-file.
Wherever relevant, I have flagged appropriate lines and mentioned the line
numbers here.
I have left detailed comments in the files themselves.
I have also added general comments, if any, over here.

## coati_function_library_V1.R

### latlon.to.utm(...)

- The argument `utm.zone` defaults to 34, where the KRR is. Pay close
      attention to this while using it in the dataset from Panama (which falls
      in 16N and 17N)! Also keep in mind the default value for the argument
      `southern_hemisphere` is `TRUE`.
- l30 flagged.

### utm.to.latlon(...)

- Same concerns as the argument above, about default argument values. Not
      a bug precisely, but a likely source of potential bugs.

### get_subgroup_data(...)

- I want to raise a concern about the NaN policy espoused here. The
      DBSCAN algorithm is run, and subgroups determined, when even one location
      is not NaN. This would be sort of fine in cases where you have an enormous
      number of data-points, but in this case, the max number of data-points at
      any time is the number of tagged coatis. This makes the grouping
      particularly susceptible to missing individuals, especially in some
      spatial conformations of the group. Depending on how much NaNs are in the
      data, I would strongly urge going for a tougher NaN policy, i.e., not
      computing subgroup identities for a point of time if _any_ location there
      is missing. If this is not feasible because of a lot of NaNs in some
      individuals, some sensitivity analyses trying out varying values of the
      DBSCAN radius might become necessary.

### visualize_network_matrix(...)

- checked

### visualize_network_matrix_trago(...)

- checked
- This seems, however, to be more or less a duplicate of the function above.
      Why not implement this as an argument in the above function instead? Much
      easier to avoid small bugs that way.

### get_proximity_data(...)

- function generates an overall-proximity matrix
- I am not entirely sure why the loop over time is necessary. You can call
      arithmetic operators on vectors pretty easily, e.g., `ds <- sqrt(xs[i,]^2
      + xs[j,]^2)`. A similar point about pointless computation is that i and
      j are both looped over all individuals. This isn't necessary because
      proximity and distances are both symmetric measures. You could therefore
      use

      ```
      for(i in 1:n_ind){
        for(j in i+1:n_ind){
            ...
        }
      }
      ```
- However, of course, the above points don't get in the way of the
      correctness of the results. This function therefore checks out.

### randomise_splits(...)

- checked
- TODO: confirm what n_sub1, n_sub2, and n_sub3 are

### dist_to_0_or_1(...)

- checked
- but could also simply be implemented as
  ```
  return((2*q) %/% 1)
  ```
  :p

### get_consistency(...)

- checked code
- but I don't see why you need the `dist_to_0_or_1` function in the first place.
  why not define your consistency simply as `1 - 2*mean(p_dyad_together)`? I did
  not find a good justification in the manuscript as well, and it seems to get
  rid of a lot of useful information (unless I'm missing something).
