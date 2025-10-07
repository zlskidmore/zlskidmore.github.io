---
layout: default
title: Contact
permalink: /contact/
---

<!-- Contact Section -->
<section class="py-5">
  <div class="container">
    <form action="https://formspree.io/f/mqaydqal" method="POST" id="contact-form" class="mx-auto" style="max-width: 600px;">
      
      <div class="mb-3">
        <label for="name" class="form-label">Name</label>
        <input type="text" name="name" id="name" class="form-control" required>
      </div>

      <div class="mb-3">
        <label for="email" class="form-label">Email</label>
        <input type="email" name="_replyto" id="email" class="form-control" required>
      </div>

      <div class="mb-3">
        <label for="message" class="form-label">Message</label>
        <textarea name="message" id="message" class="form-control" rows="5" required></textarea>
      </div>

      <button type="submit" class="btn btn-primary w-100">Send Message</button>

      <!-- Success/Failure Message Placeholder -->
      <div id="form-status" class="mt-3 text-center"></div>
    </form>
  </div>
</section>
