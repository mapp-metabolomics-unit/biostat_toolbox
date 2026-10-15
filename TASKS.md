# Tasks

This is the user-owned task inbox. Append tasks under **Queue**; the agent must not edit this file or `ANSWERS.md`. The agent records progress, completion, and questions in `AGENT_LOG.md`, using task IDs. Answer questions in `ANSWERS.md` using the question IDs from the log.

- Give each task a permanent ID (`T003`, `T004`, …); never renumber. Subtasks use `T003.1`, `T003.2`, etc.
- Put searchable labels in brackets, such as `[dashboard] [bug]`. Labels are optional.
- No checkboxes or status edits here: consult `AGENT_LOG.md` for status. When blocked, the agent should log a numbered question, work on another runnable task, and re-read `ANSWERS.md` before the next pass. Re-read the queue before finishing. A finished agent still needs a new message to resume.

Copy this shape when adding a task (omit subtasks or acceptance criteria if unnecessary):

```text
### T003 [area] [type] Short title
What should change, and any relevant context.

- T003.1 First distinct part
- T003.2 Second distinct part

Done when: Observable result or acceptance criteria.
```

## Completed before this queue

### T001 [dashboard] Tables
For each table displayed in the dashboard:
- T001.1 It can be downloaded.
- T001.2 Each column can be sorted and filtered.
- T001.3 Displayed columns can be selected from a dropdown menu.
- T001.4 When there is a numerical value in a field. How can the user filter for a range or above a threshold ? e.g. all Saved InChIKey2D agreement abbove 2 ? Implement a solution

### T002 [dashboard] Volcano plots
- T002.1 Center the plot; the log2 fold-change line should always be centered.

## Queue

<!-- Append new tasks below this line. Start with T003. -->

### T003 [dashboard] UI

- T003.1 When hovering a feature related dataplot, display relevant metadata of the feature. Typically feature id, mz, rt. If a chemical structure exists, display it in the hover box. 
- T003.2 Add Settings menu on the left
    - T003.2.1 In this menu add possibility ofor the user to select which fields should be display in the hover menu
- T003.3 Rename "Saved source structures · not confirmed identities 2D images load from the Natural Products Cheminformatics API, as in the Shiny viewer. This sends each displayed SMILES and your network address to that service. Distinct saved SMILES are grouped by exact text; no identity is inferred." o something simpler. E.g Chemical depiction + a link to the API  No need to refer to Shiny.
        Done when: User can browse and hover a volcano plot and view chemical structure without clicking. He can select which source can be used for the SMILES display
- T003.4 overall make sure to homogeneize sources of structures SMILES. We have 3 possibilities GNOPS, ISDB and SIRIUS (met_annot_enhancer reports ISDB structure) it is thus redudnat. Make sure this is 
- T003.5 Make the foldable feature details bars sliglty wider.
- T003.6 Make the foldable feature exploration panel detachable. So We can move it to an external screen and keep it synced wihile exploring data. If not possible report in the log.
- T003.7 When I select specific setting for the feature boxplot I expect these settings to be saved and persistent when selecting another feature. Atm, it is not the case. Fix.



### T004 [dashboard] Filters

- T004.1 The filters on the left bar such as Colour/ group by Sample role Groups to show. Should dissapear since they are redundant with the filetring option in the PCA, PCoA Volcano etc. plots.
- T004.2 When the filetrs are not redundant for example under Descriptive sample intensities, Then the grouping option should appear there (and not in the left side bar.)

### T005 [dashboard] Feature boxplot

- T005.1 Add facet options
- T005.2 Add color picker 
- T005.2.1 The color picker for groups should be externalized and made general for the dashboard. The user should not nbe able to fix different color for the same group according to the plots- 
- T005.3 Hide the Saved intensity stage function s"Group by Sample role Facet by Groups to show Ctrl/Cmd-click to show multiple groups. These controls affect only this descriptive plot. Group colours" So they are accesible only when unfolding a display option menu. We want to Keep the UX as smooth as possible

### T006 [dashboard] PCA

- T006.1 In theFeature loadings plot have a color code for annotated versus not tannotated features (annotated means at least one SMILES in GNPS, ISDB or SIrius)
- T006.2 Sync the display of Feature loadings and the PCA. In 3D mode when mooving one plot loading or pca the other one should move in a synced manner. If too complex or hacky do not implement. Justr a nice to have.


### T007 [dashboard] Metadata

- T007.1 Add a metadata viewer section under the Overview tab.

### T008 [stats_workflow] rerun

- T008.1 Rreun a workflow with a PLS-DA (a supervised analysis) of WT_O versus xdhB2.
- T008.2 Once T010.1 is done. Rerun a stats workflow with the new ropls functions.

### T009 [dashboard] features

- T009.1 overall the dashboard for all plots displaying multiple features (e.g. loading of PCA, Volcano) add a foldable Display options menu allowing to mapp a given column with a color code to the features. For ex. on the loading of the PCA. I can mapp a color for each CANOPUS NPC pathway  associated to the features.
- T009.2 when clicking a feature in any plot it should Open a foldable version. of the Feature explorer panel so as to allow the user to navigate on the same screen on bothe for example the Volcano plot and the the feature barplots + metadata
- T009.3 An annotated feature should only display SMILES  when IK of sources are different,. That is if ISDB, GNPS and SIRIUS agree we should display one unique SMILES.
- T009.4 Implement a 3 Star system to be displayed near a structure when it has 1, 2 or 3 tools agreeing. One of the Stars (corresponding to ISDB) could be color coded ranging from pale yellow to bright green when the taxonomical score increases. In other words an annotation with 3 stars and one of them bright green  would mean in this case consensius annotation reported in Ensifer meilotus.
- T009.5 The feature filters hereunder Color feature by should also offer to color feature by number of consensus annotation (o, 1, 2 or 3) at m the candidate record source is unclear.

### T010 [stats_workflow] logic

- T010.1 We use already published and verified implementation of the unsupervised and supervised statistics. See https://bioconductor.org/packages/release/bioc/html/ropls.html 
- T010.2 For the PCoA, if applicable, I believe we should also have a feature loadings plot
- T010.3 Just as for PCA, for the supervised statistics I would expect o have loadings plots and even S plots or VIP plots if possible see https://bioconductor.org/packages/3.23/bioc/vignettes/ropls/inst/doc/ropls-vignette.html#44_Orthogonal_partial_least_squares:_OPLS_and_OPLS-DA
- T010.4 Create a new tab that guides the user for the creation of a recipe for the stats. Example It should allow me from the current metadata table to say that I want to build a new group gathering WT_N and WT_O versus all xdhA (1 and 2) versus alll xdhB. This should feed from the general metadata file. Then I could select what type of filters I want, what type of stats to be run etc. etc-. Ok to have it in an Advanced tab..
- T010.5 Not happy with the way we store the outputs in /home/allardpm/git_repos/mapp-metabolomics-unit/didier-reinhardt-group/docs/mapp_project_00007/mapp_batch_00275/results/stats_v2. First why do we have /home/allardpm/git_repos/mapp-metabolomics-unit/didier-reinhardt-group/docs/mapp_project_00007/mapp_batch_00275/results/stats_v2_inputs and duplicates of the metadata, quant table and all. There should be one unique source for these files and they should not be altered asfaict. Then under /home/allardpm/git_repos/mapp-metabolomics-unit/didier-reinhardt-group/docs/mapp_project_00007/mapp_batch_00275/results/stats_v2 we no longer have the /home/allardpm/git_repos/mapp-metabolomics-unit/didier-reinhardt-group/docs/mapp_project_00007/mapp_batch_00275/results/stats/params_log.tsv a mechanisms as such (and the list of jobs plus a short human readable description should be something that can be browsed fromn the dashboard.)
- T010.5 In /home/allardpm/git_repos/mapp-metabolomics-unit/didier-reinhardt-group/docs/mapp_project_00007/mapp_batch_00275/recipes/merged.recipe.yaml how do you expect me to set analysis_group as target group ? Wher is that setup ?

### T011 [dashboard] chemical class exploration

- T011.1 based on the same logic we had established for the /home/allardpm/git_repos/mapp-metabolomics-unit/biostat_toolbox/src/generate_data_explorer.R add a tab Allowing to drill down features starting from their NPC pathways then super class then cklass then feature. The groups to be compared should be customizable 

### T012 UX

- T012.1 When running commands such as rig run -r 4.6.1 -f run.R -- run --dataset
../../didier-reinhardt-group/docs/mapp_project_00007/mapp_batch_00275/recipes/merged.dataset.yaml --recipe
../../didier-reinhardt-group/docs/mapp_project_00007/mapp_batch_00275/recipes/merged.recipe.yaml
The users would like more verbosity and ideally a progress bar for long steps.
