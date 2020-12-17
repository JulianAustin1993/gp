name := "gp"

version := "0.1"

scalaVersion := "2.13.4"

libraryDependencies += "org.scalatest" %% "scalatest" % "3.0.8" % Test

libraryDependencies ++= Seq(
  // Last stable release
  "org.scalanlp" %% "breeze" % "1.1",

  // Native libraries are not included by default. add this if you want them
  // Native libraries greatly improve performance, but increase jar sizes.
  // It also packages various blas implementations, which have licenses that may or may not
  // be compatible with the Apache License. No GPL code, as best I know.
  "org.scalanlp" %% "breeze-natives" % "1.1",

  // The visualization library is distributed separately as well.
  // It depends on LGPL code
  "org.scalanlp" %% "breeze-viz" % "1.1"
)

libraryDependencies += "org.typelevel" %% "spire" % "0.17.0-M1"


