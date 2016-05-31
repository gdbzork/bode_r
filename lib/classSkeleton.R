# Generic skeleton of a class with two subclasses

setGeneric("zork",def=function(obj) standardGeneric("zork"))
setClass("Base",
         representation(field1="character", field2="numeric"),
         prototype(field1="Zork",field2=4)
)

setValidity("Base",function(object) nchar(object@field1) == object@field2)
setMethod("show","Base",function(object)
                          print(paste("Base:",
                                      object@field1,
                                      "(len=",
                                      object@field2,
                                      ")")))
setMethod("zork","Base",function(obj) { print(obj@field2) })
setMethod("initialize","Base",function(.Object, x="zork") {
                                .Object@field1 <- x
                                .Object@field2 <- nchar(x)
                                .Object })

setClass("Sub1",
         contains=c("Base"),
         representation=representation(field3="numeric"),
         prototype(field3=99)
)
setMethod("show","Sub1",function(object)
                          print(paste("Sub1:",
                                      object@field1,
                                      "(len=",
                                      object@field2,
                                      ", other=",
                                      object@field3,
                                      ")")))

setMethod("initialize","Sub1",function(.Object) {
                                .Object@field3 = - .Object@field2
                                .Object})


setClass("Sub2",
         contains=c("Base"),
         representation=representation(field3="numeric"),
         prototype(field3=-99)
)

setMethod("show","Sub2",function(object)
                          print(
                            paste("Sub2:",
                                  object@field1,
                                  "(len=",
                                  object@field2,
                                  ", other=",
                                  object@field3,
                                  ")")))

setMethod("initialize","Sub2",function(.Object) {
                                .Object@field3 = .Object@field2^2
                                .Object})

