## Personal Website — Zachary L. Skidmore
Development repository for the personal and professional website of Zachary L. Skidmore.
This site is built with [Jekyll](https://jekyllrb.com) - 4.4.1, using Ruby - 3.4.4 and Bundler - 2.7.2.

### Development Setup
To develop locally, the recommended workflow is:

#### 1. Install Ruby via rbenv
Follow the rbenv [installation guide](https://github.com/rbenv/rbenv).

Once installed, run in the site directory:
```
rbenv install 3.4.4
rbenv local 3.4.4
```

#### 2. Install Bundler and Jekyll
```
gem install bundler jekyll
bundle install
```

#### 3. Run the development server
```
bundle exec jekyll serve --livereload
````
*site will be available on localhost:4000*

### Feedback
Questions, comments, or suggestions about the site or its content?
Please open a [GitHub Issue](https://github.com/zlskidmore/zlskidmore.github.io/issues) in this repository.

