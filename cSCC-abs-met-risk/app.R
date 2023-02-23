#---------------------------------------------------
# Aim: Deploy R shiny app so that users can use the model with a patient
# Author: B. Rentroia Pacheco
# Input: Patient Characteristics
# Output: Prediction of metastatic risk
#---------------------------------------------------


# Load packages ----
library(shiny)

# Source helpers ----
# function to compute the linear predictor of the model, based on the patient and tumor characteristics provided by the user.
predict_lp = function(input){
  pat_characteristics = unlist(input)
  #print(pat_characteristics)
  # Model
  #coefs_model =c(0.25467,0.51388,0.574366,-0.72485,0.51388,0.540523,0.327348,1.401381,1.296865,0.801249)
  coefs_model =c(0.25,0.51,0.57,-0.72,0.52,0.54,0.33,1.4,1.3,0.8)
  names(coefs_model)=c("age","sex","n_cscc","t_location_sneck","t_location_face","t_diam","tinv_sc","tinv_beysc","differentiation","pni_or_lvi")
  
  # Adapt input
  pat_characteristics["t_diam"]= ifelse(as.numeric(pat_characteristics["t_diam"])>4,4,pat_characteristics["t_diam"])
  pat_characteristics["n_cscc"] = ifelse(as.numeric(pat_characteristics["n_cscc"])>4,4,as.numeric(pat_characteristics["n_cscc"]))
  pat_characteristics["t_location_sneck"] = ifelse(pat_characteristics["t_location"]=="Scalp/Neck",1,0)
  pat_characteristics["t_location_face"] = ifelse(pat_characteristics["t_location"]=="Face",1,0)
  pat_characteristics["tinv_sc"] = ifelse(pat_characteristics["t_involvement"]=="SubFat",1,0)
  pat_characteristics["tinv_beysc"] = ifelse(pat_characteristics["t_involvement"]=="BeySubFat",1,0)
  # Remove mean values of continuous variables
  pat_characteristics["age"]= (as.numeric(pat_characteristics["age"]))/10 - 7.47
  pat_characteristics["n_cscc"]= as.numeric(pat_characteristics["n_cscc"])-0.29 
  pat_characteristics["t_diam"]= as.numeric(pat_characteristics["t_diam"]) - 1.20 
  
  pat_characteristics=as.numeric(pat_characteristics[names(coefs_model)])
  
  
  # compute linear predictor
  lp = sum(as.numeric(t(coefs_model)%*%pat_characteristics))
  
  return(lp)
}

#Function to compute the metastatic risk at timepoint tm, given the patient and tumor characteristics provided by the user.
predict_met_risk = function(tm,input){
  
  lp = predict_lp(input)
  # Baseline survival at different time points
  baseline_survival = c(0.987,0.976,0.973)
  names(baseline_survival)=as.character(c(1,3,5))
  tm = as.character(tm)
  
  mean_lp = 1.86
  
  # compute metastatic risk
  met_risk = round(100*(1-baseline_survival[tm]^exp(lp-mean_lp)),1)
  
  return(met_risk)
}

# User interface ----
ui <- fluidPage(
  tags$head(
    tags$style(
      HTML(".shiny-notification {
             position:fixed;
             top: calc(50%);
             left: calc(50%);
             }
             "
      )
    )
  ),
  
  # App title ----
  titlePanel(
    fluidRow(
      column(10, "Prediction of metastatic risk in patients with cutaneous squamous cell carcinoma (cSCC)",    
             h4("This web-based calculator has been developed by the Skin Cancer Research Group of the Department of Dermatology at the Erasmus MC Cancer Institute in Rotterdam, and validated in an external cohort of cSCC patients, as described in this publication [Add link]."),
             h4("This model has only been developed and validated in cSCC and not in mucosal or genital SCC."),
      ),       column(2, img(height = 80, width = 170,src='EMC_Cinst_logo.png', align = "right"))),
    windowTitle = "EMC cSCC absolute risk model"),
  # Sidebar layout with input and output definitions ----
  # No sidebar is needed.
    
    
  # Main panel for displaying outputs ----
   mainPanel(
      fluidRow(
        column(12,h4(strong("Patient characteristics")),
               column(3,numericInput("age", label = "Age, in years", value = "-",min=18,max=120),
                helpText("Age, in years. Patients should be adults (18 years or older)")),
               column(3,selectInput("sex", label = "Sex", 
                           choices = list(" "=NA,"Female" = 0, "Male" = 1), 
                           selected = NA)),
               column(3,selectInput("n_cscc", label = "Number of prior cSCCs", 
                           choices = list(" "=NA,"0" = 0, "1" = 1, "2" = 2, "3" = 3, "4" = 4,"More than 4"=5), 
                           selected = NA)),
               column(3,selectInput("t_location", label = "Tumor location", 
                           choices = list(" "=NA,"Trunk and extremities" = "Trunk/Ext", "Scalp or neck" = "Scalp/Neck", "Face (including ear)" = "Face"), 
                           selected = NA))
        ),
        column(12,h4(strong("Tumor characteristics")),
               column(3,numericInput("t_diam", label = "Tumor diameter, in cm", value = NA,min=0,max=40,step = 1),
                      helpText("Macroscopic tumor diameter, as measured by a pathologist, in centimeters. Decimals are allowed e.g. 2.2")),
               column(3,selectInput("t_involvement", label = "Tissue involvement", 
                           choices = list(" "=NA,"Dermis" = "Dermis", "Subcutaneous fat" = "SubFat", "Beyond subcutaneous fat" = "BeySubFat"), 
                           selected = NA),
               helpText("Deepest layer of tissue involvement")),
               column(3,selectInput("differentiation", label = "Differentiation Grade", 
                           choices = list(" "=NA,"Good/moderate" = 0, "Poor/undifferentiated" = 1), 
                           selected = NA),
                  helpText("Differentiation grade according to the adjusted Broder classification system: good/moderate differentiation 0-75% undifferentiated cells, poor/undifferentiated >75% undifferentiated cells.")),
               column(3,selectInput("pni_or_lvi", label = "Perineural or lymphovascular invasion", 
                           choices = list(" "=NA,"Absent" = 0, "Present" = 1), 
                           selected = NA),
                  helpText("Presence of perineural invasion (>=0.1mm) or lymphovascular invasion of any size"))),
      width=12),
      actionButton("action", label = "Predict",width="100%",class = "btn-primary btn-lg"),width=23),
  br(),
  br(),
  br(),
  column(
    h4("A patient with these characteristics has a probability of :"),
    h4(uiOutput("myList")),
    width=12),
    br(),
    br(),
    br(),
    br(),
    br(),
    p(strong("DISCLAIMER:"),"Currently, this calculator can only be used for research purposes. It cannot be used as a medical device. Erasmus MC has no liability for any medical decision taken based on the information provided by the calculator."),
)

# Server:
server <- function(input, output) {
  
  # Global variables to keep track of the changes in the session
  # Disclaimer variable is meant to keep track of whether the user has agreed with the disclaimer, so that the s/he only have to agree with it once.
  disclaimer <<- reactiveValues()
  disclaimer$value= FALSE
  # Button predict variable is meant to only change the risk probabilities if the user clicks on the predict button.
  # #n will be updated each time the input variables change, but the risk probability will only be displayed if 
  button_predict <<-list()
  button_predict$n <<- 0
  
  #This will update the metastatic risk probabilities when the user has accepted the disclaimer and pressed the submit button.
  # We know the submit button has been pressed because in that case, button_predict$n = input$action.
  met_risk_interactive = reactive({
    if (disclaimer$value & (button_predict$n ==input$action)){
      button_predict$n <<- button_predict$n+1
      
      #Check numerical inputs are correct:
      if(input$age<18 | input$age>120){
        showModal(modalDialog(
          title = p("Age should be between 18 and 120 years old."),
          easyClose = TRUE,
        ))
        c(" "," "," ")
      }else if(input$t_diam<=0){
        showModal(modalDialog(
          title = p("Tumor diameter should be larger than 0 cm."),
          easyClose = TRUE,
        ))
        c(" "," "," ")
      }else{
        sapply(c(1,3,5),predict_met_risk,reactiveValuesToList(input))
      }
      
    }else{
      button_predict$n <<- button_predict$n+1
      c(" "," "," ")
    }})
  
  #Display metastatic risk.
  output$myList <- renderUI(HTML(markdown::renderMarkdown(text = paste(paste0("- **",met_risk_interactive() ,"** % of developing a metastasis within ",c(1,3,5),c(" year",rep(" years",2)), "\n"), collapse = "")
  )))
  
  #Actions after predict button has been clicked on
  observeEvent(input$action, {
    #All variables should be filled:
    if (sum(reactiveValuesToList(input)%in%c("NA",NA))==0){
      #Show the disclaimer (only once))
      if(disclaimer$value==FALSE){
        showModal(modalDialog(
          title = p(strong("DISCLAIMER:"),"Currently, this calculator can only be used for research purposes. It cannot be used as a medical device. Erasmus MC has no liability for any medical decision taken based on the information provided by the calculator."),
          easyClose = TRUE,
          footer = modalButton("I agree")
        ))
        disclaimer$value=TRUE
      }
      button_predict$n <<- input$action #This will enable the risk probability to be displayed.
      met_risk_interactive()
    }else{
      showModal(modalDialog(
        title = p("Please fill in all variables"),
        easyClose = TRUE,
      ))
    }
  })
  
  
}

# Run the app
shinyApp(ui, server)

# Note: This is the code to upload the app in shinyapps.io:
#rsconnect::deployApp('C:/Users/b.rentroia/Documents/CP_refined Repository/cSCC-abs-met-risk/cSCC-abs-met-risk',account="emc-dermatology")
