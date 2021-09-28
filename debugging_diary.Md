# Helen's notes to self
## 09/27/2021
- Having trouble with the f strings of the regulation term
- The regulation term will NOT give me the equation for the new transcription rate alone, only as part of the gRNA equation.
    - Jake says having the Hill equaiton only embedded is fine, solved
- How do I get rid of the original transcription rate so that I only have the new transcription rate in the formula?
    - Jake said: Don't put the base rate in the equation, solved
- Prettied up the equations
    - Got the frac to be actually split level had to do "\\frac{{a}}{{b}}"
        - Or use triple brackets if I'm calling a variable
## 09/26/2021
- I got rid of the 2nd degradation term by removing the explicit "degradation" interaction in the model creation (It didn't look like cas had one)
- I wasn't able to get into the "regulation_term" function that Jake pre-made for me, just nothing was satisfying the conditions
    - The only thing that had a "regulation" list WASN't a SBO_PRODUCT
- I was able to get into that function by changing my make model code, I hade to have TF acting on the sgRNA rna product, NOT the sgRNA DNA that Jake was returning
    - Not sure if I can move up or down a level, so for now I have the function returning both the RNA and the DNA but am only using the RNA downstream
- I am re-defining i_type in that function (don't want to be redundant in the variables I am passing) and then if statements the possible ones I could see
- Not entirely sure what a regulation term is though
- Not sure why the regulation term is acting on V
## Pre 09/26/2021
- I was generating a model doc
- I added the name to the definitions of the LaTeX so equations were being generated
- But the equations don't have the interaction between the TF and the RNA
- And there are two degradation terms for the TF